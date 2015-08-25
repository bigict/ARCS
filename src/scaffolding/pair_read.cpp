#include "pair_read.h"
#include "kseq.h"

#include <fstream>
#include <numeric>
#include <set>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffording.pair_read"));

class PairReadReader {
public:
    PairReadReader(std::istream& stream1, std::istream& stream2) : _reader1(stream1), _reader2(stream2) {
    }

    bool read(PairRead& read) {
        DNASeq read1, read2;
        if (_reader1.read(read1) && _reader2.read(read2)) {
            read.set1 = read1.seq;
            read.set2 = read2.seq;
            return true;
        }
        return false;
    }
private:
    DNASeqReader _reader1;
    DNASeqReader _reader2;
};

typedef boost::accumulators::accumulator_set< size_t, boost::accumulators::stats< boost::accumulators::tag::mean, boost::accumulators::tag::moment< 2 > > > Accumulator;
struct AccuWrapper {
    AccuWrapper(Accumulator& acc, size_t min_threshold) : _acc(acc), _min_threshold(min_threshold) {
    }
    void operator()(size_t item) {
        if (item > _min_threshold) {
            _acc(item);
        }
    }
private:
    Accumulator& _acc;
    size_t _min_threshold;
};

void InsertSizeEstimater::estimate(size_t* mean, double* variance) const {
    InsertSizeDistr insert_size_distr;

    BOOST_FOREACH(const PairRead& pair_read, _pair_reads) {
        if (insert_size_distr.size() >= 1000) {
            break;
        }
        if (pair_read.set1.size() < _K || pair_read.set2.size() < _K) {
            continue;
        }

        distribution(pair_read.set1, make_complement_dna(pair_read.set2), insert_size_distr);
        distribution(pair_read.set2, make_complement_dna(pair_read.set1), insert_size_distr);
    }

    Accumulator acc;
    std::for_each(insert_size_distr.begin(), insert_size_distr.end(), AccuWrapper(acc, 10));

    if (mean != NULL) {
        *mean = boost::accumulators::mean(acc);
    }
    if (variance != NULL) {
        *variance = boost::accumulators::moment< 2 >(acc) - std::pow(boost::accumulators::mean(acc), 2);
    }
}

void InsertSizeEstimater::distribution(const std::string& read1, const std::string& read2, InsertSizeDistr& insert_size_distr) const {
    LOG4CXX_TRACE(logger, boost::format("pair read %s %s") % read1 % read2);

    Kmer kmer1, kmer2;
    std::pair< size_t, long > left, right;
    for (size_t i = 0, j = _K; j < read1.size() && j < read2.size(); ++i, ++j) {
        kmer1 = read1.substr(i, _K);
        kmer2 = read2.substr(i, _K);

        left = _kmer_tbl.findPos(kmer1);
        right = _kmer_tbl.findPos(kmer2);
        LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre1=%s %d %d") % kmer1 % left.first % left.second);
        LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre2=%s %d %d") % kmer2 % right.first % right.second);

        if (left.first != -1  && right.first != -1 && left.first == right.first) {
            if (right.second - left.second > _insert_size / 4 
                    && right.second - left.second < 2 * _insert_size) 
                insert_size_distr.push_back(right.second - left.second);
        }
    }
}

PairReadSet::PairReadSet(std::istream& stream1, std::istream& stream2, size_t K, size_t INSERT_SIZE) : _k(K) , INSERT_SIZE(INSERT_SIZE) {
    DELTA = 0;
    PAIR_KMER_NUM = 0;
    
    init(stream1, stream2);
}

PairReadSet::PairReadSet(const std::string& file1, const std::string& file2, size_t K, size_t INSERT_SIZE) : _k(K) , INSERT_SIZE(INSERT_SIZE) {
    DELTA = 0;
    PAIR_KMER_NUM = 0;
    
    std::ifstream stream1(file1.c_str()), stream2(file2.c_str());
    init(stream1, stream2);
}

void PairReadSet::init(std::istream& stream1, std::istream& stream2) {
   LOG4CXX_INFO(logger, boost::format("read pair read begin"));

   PairReadReader reader(stream1, stream2);
   PairRead pair_read;
   while (reader.read(pair_read)) {
        _pair_reads.push_back(pair_read);
   }

   LOG4CXX_INFO(logger, boost::format("read pair read end"));
}

void PairReadSet::buildConnectGraph(Graph& g, KmerTable& tbl, const ContigSet& contigset, const std::vector<Component>& components) {
    LOG4CXX_INFO(logger, boost::format("in build graph: K=%lld") % _k);

    std::string read1, read2;
    BOOST_FOREACH(const PairRead& pair_read, _pair_reads) {
        const std::string& read1 = pair_read.set1;
        const std::string& read2 = pair_read.set2;
        if (read1.size() < _k || read2.size() < _k) {
            continue;
        }
        findLink(read1, make_complement_dna(read2), g, tbl, contigset, components);
        findLink(read2, make_complement_dna(read1), g, tbl, contigset, components);
    }
}

void PairReadSet::estimateInsertSize(const KmerTable &tbl) {
    InsertSizeEstimater estimater(_k, INSERT_SIZE, _pair_reads, tbl);
    double variance = 0;
    estimater.estimate(&INSERT_SIZE, &variance);
    DELTA = std::sqrt(variance);

    LOG4CXX_INFO(logger, boost::format("INSERT=%d DELTA=%d") % INSERT_SIZE % DELTA);
    return;

    std::vector< int > insert_len_dis;

    BOOST_FOREACH(const PairRead& pair_read, _pair_reads) {
        if (insert_len_dis.size() >= 1000) {
            break;
        }
        if (pair_read.set1.size() < _k || pair_read.set2.size() < _k) {
            continue;
        }

        estimateOnePR(pair_read.set1, make_complement_dna(pair_read.set2), tbl, insert_len_dis);
        estimateOnePR(pair_read.set2, make_complement_dna(pair_read.set1), tbl, insert_len_dis);
    }

    size_t sum = 0;
    size_t index = 0;
    BOOST_FOREACH(const int insert_num, insert_len_dis) {  
        if (insert_num > 10) {
            sum += insert_num;
            index ++;
        }
    }
    if (index > 0) {
        INSERT_SIZE = double(sum) / index;
        sum = 0;
        BOOST_FOREACH(const int insert_num, insert_len_dis) {
            if (insert_num > 10)
                sum += (insert_num - INSERT_SIZE)*(insert_num - INSERT_SIZE);
        }
        DELTA =  std::sqrt(double(sum) / index);
    }

    LOG4CXX_INFO(logger, boost::format("INSERT=%d DELTA=%d") % INSERT_SIZE % DELTA);
}

void PairReadSet::estimateOnePR(const std::string& read1, const std::string& read2, const KmerTable& kmer_tbl, std::vector<int>& insert_length_list) {
    LOG4CXX_TRACE(logger, boost::format("pair read %s %s") % read1 % read2);

    Kmer kmer1, kmer2;
    std::pair<size_t, long> left,right;
    for (int i=0; i<read1.size()-_k+1 && i<read2.size()-_k+1;i++) {
        kmer1 = read1.substr(i, _k);
        kmer2 = read2.substr(i, _k);

        left = kmer_tbl.findPos(kmer1);
        right = kmer_tbl.findPos(kmer2);
        LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre1=%s %d %d") % kmer1 % left.first % left.second);
        LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre2=%s %d %d") % kmer2 % right.first % right.second);

        if (left.first != -1  && right.first != -1 && left.first == right.first) {
            if (right.second - left.second > INSERT_SIZE / 4 
                    && right.second - left.second < 2*INSERT_SIZE) 
                insert_length_list.push_back(right.second - left.second);
        }
    }

}

struct PairCmp{
    bool operator()(const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b) {
        if (a.first != b.first) {
            return a.first < b.first;
        }
        return a.second < b.second;
    }
};

void PairReadSet::findLink(const std::string& read1, const std::string& read2, Graph& graph, const KmerTable& tbl, const ContigSet& contigset, const std::vector<Component>& components) {
    std::set< std::pair< size_t, size_t >, PairCmp > find_component;
    for (int i = 0; i < read1.size() - _k + 1 && i < read2.size() - _k + 1; ++i) {
        Kmer kmer1 = read1.substr(i, _k);
        Kmer kmer2 = read2.substr(i, _k);

        std::pair<size_t, long> left = tbl.findPos(kmer1);
        std::pair<size_t, long> right = tbl.findPos(kmer2);
        
        if (left.first != right.first && left.first != -1 && right.first != -1) {
            PAIR_KMER_NUM ++;

            LOG4CXX_TRACE(logger, boost::format("pair kmre1=%s %d %d") % kmer1 % left.first % left.second);
            LOG4CXX_TRACE(logger, boost::format("pair kmre2=%s %d %d") % kmer2 % right.first % right.second);
        
            long len_tmp = INSERT_SIZE + left.second - right.second;
            long left_len = Component::length(_k, contigset, components[left.first]);
            long overlap = left_len - len_tmp - _k + 1;
            if (overlap > 0) {
                if (overlap > INSERT_SIZE)
                    continue;
                len_tmp = left_len - _k + 1;
            }
            if (find_component.find(std::make_pair(left.first, right.first)) != find_component.end()) {
                graph.addEdge(left.first, right.first, len_tmp, false);
            } else {
                find_component.insert( std::make_pair(left.first, right.first));
                graph.addEdge(left.first, right.first, len_tmp, true);
            }
        }
    }
}

std::ostream& operator<<(std::ostream& os, const PairReadSet& p_r) {
    for (int i=0; i<p_r._pair_reads.size(); i++) {
        os << p_r._pair_reads[i].set1 <<"\t"<<p_r._pair_reads[i].set2 << std::endl; 
    }
    return os;
}

bool ReadPairReads(std::istream& stream1, std::istream& stream2, PairReadList& pair_reads) {
   LOG4CXX_INFO(logger, boost::format("read pair read from stream begin"));

   PairReadReader reader(stream1, stream2);
   PairRead pair_read;
   while (reader.read(pair_read)) {
        pair_reads.push_back(pair_read);
   }

   LOG4CXX_INFO(logger, boost::format("read pair read from stream end"));
}

bool ReadPairReads(const std::string& file1, const std::string& file2, PairReadList& pair_reads) {
    std::ifstream stream1(file1.c_str()), stream2(file2.c_str());
    return ReadPairReads(stream1, stream2, pair_reads);
}

