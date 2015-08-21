#include "pair_read.h"
#include "kseq.h"

#include <fstream>
#include <numeric>
#include <set>

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

void PairReadSet::buildConnectGraph(Graph& g, KmerTable& tbl, const std::vector<Component>& components) {
    std::string read1, read2;
    std::vector<int> insert_len_dis;
    LOG4CXX_INFO(logger, boost::format("in build graph: K=%lld") % _k);
    BOOST_FOREACH(const PairRead& pair_read, _pair_reads) {
        read1 = pair_read.set1;
        read2 = pair_read.set2;
        if (read1.size() < _k || read2.size() < _k)
            continue;
        std::string read1_inv = make_complement_dna(read1);
        std::string read2_inv = make_complement_dna(read2);
        findLink(read1, read2_inv, g, tbl, components);
        findLink(read2, read1_inv, g, tbl, components);
    }
}

struct InsertLengthPlus {
    InsertLengthPlus(size_t threshold) : _threshold(threshold) {
    }
    size_t operator()(size_t l, int r) const {
        return r > _threshold ? l + r : l;
    }
private:
    size_t _threshold;
};

void PairReadSet::estimateInsertSize(const KmerTable &tbl) {
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

void PairReadSet::findLink(const std::string& read1, const std::string& read2, Graph& graph, const KmerTable& tbl, const std::vector<Component>& components) {
    std::set< std::pair<size_t, size_t>, PairCmp > find_component;
    for (int i=0; i<read1.size()-_k+1 && i<read2.size()-_k+1; i++) {
        Kmer kmer1 = read1.substr(i, _k);
        Kmer kmer2 = read2.substr(i, _k);

        std::pair<size_t, long> left = tbl.findPos(kmer1);
        std::pair<size_t, long> right = tbl.findPos(kmer2);
        
        if (left.first != right.first && left.first != -1 && right.first != -1) {
            PAIR_KMER_NUM ++;
            LOG4CXX_TRACE(logger, boost::format("pair kmre1=%s %d %d") % kmer1 % left.first % left.second);
            LOG4CXX_TRACE(logger, boost::format("pair kmre2=%s %d %d") % kmer2 % right.first % right.second);
        
            long len_tmp = INSERT_SIZE + left.second - right.second;
            long left_len = components[left.first].length();
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
