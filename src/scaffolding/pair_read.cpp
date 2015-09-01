#include "pair_read.h"
#include "graph.h"
#include "kmer.h"
#include "kmer_tbl.h"
#include "kseq.h"

#include <fstream>
#include <numeric>
#include <set>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffolding.pair_read"));

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

bool ReadPairReads(std::istream& stream1, std::istream& stream2, PairReadList& pair_reads) {
    PairReadReader reader(stream1, stream2);
    PairRead pair_read;
    while (reader.read(pair_read)) {
        pair_reads.push_back(pair_read);
    }
    return true;
}

bool ReadPairReads(const std::string& file1, const std::string& file2, PairReadList& pair_reads) {
    std::ifstream stream1(file1.c_str()), stream2(file2.c_str());
    return ReadPairReads(stream1, stream2, pair_reads);
}

typedef boost::accumulators::accumulator_set< int, boost::accumulators::stats< boost::accumulators::tag::count, boost::accumulators::tag::mean, boost::accumulators::tag::moment< 2 > > > Accumulator;
struct Statistics {
    Statistics(int min_threshold, Accumulator& acc) : _min_threshold(min_threshold), _acc(acc) {
    }
    void operator()(int& insert_size) const {
        if (insert_size >= _min_threshold) {
            _acc(insert_size);
        }
    }
private:
    Accumulator& _acc;
	int _min_threshold;
};

void InsertSizeEstimater::estimate(size_t* insert_size, double* delta) {
    LOG4CXX_DEBUG(logger, boost::format("pair_reads=[%d]") % _pair_reads.size());

    InsertSizeDistr insert_size_distr;
    BOOST_FOREACH(const PairRead& pair_read, _pair_reads) {
        if (insert_size_distr.size() >= 1000) {//trick
            break;
        }
        if (pair_read.set1.length() < _K || pair_read.set2.length() < _K) {
            continue;
        }
        estimateOnePR(pair_read.set1, make_complement_dna(pair_read.set2), insert_size_distr);
        estimateOnePR(pair_read.set2, make_complement_dna(pair_read.set1), insert_size_distr);
    }

    Accumulator acc;
    std::for_each(insert_size_distr.begin(), insert_size_distr.end(), Statistics(11, acc));
    if (insert_size != NULL) {
        *insert_size = boost::accumulators::mean(acc);
    }
    if (delta != NULL) {// delta = E(xi-mean)^2
        size_t mean = (size_t)boost::accumulators::mean(acc);
        *delta = std::sqrt(
                boost::accumulators::moment< 2 >(acc) - 2 * boost::accumulators::mean(acc) * mean + std::pow(mean, 2)
                );
    }
}

void InsertSizeEstimater::estimateOnePR(const std::string& read1, const std::string& read2, InsertSizeDistr& insert_size_distr) {
    LOG4CXX_TRACE(logger, boost::format("InsertSizeEstimater::estimateOnePR %s %s") % read1 % read2);

    for (size_t i = 0, j = _K; j <= read1.size() && j <= read2.size(); ++i,++j) {
        Kmer kmer1 = read1.substr(i, _K), kmer2 = read2.substr(i, _K);
        KmerList::const_iterator left = _hash_tbl.find(kmer1), right = _hash_tbl.find(kmer2);

        if (left != _hash_tbl.end() && right != _hash_tbl.end() && left->second.first == right->second.first) {
            LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre1=%s %d %d") % kmer1 % left->second.first % left->second.second);
            LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre2=%s %d %d") % kmer2 % right->second.first % right->second.second);

            if (right->second.second - left->second.second > _insert_size / 4 
                    && right->second.second - left->second.second < 2 * _insert_size) {
                insert_size_distr.push_back(right->second.second - left->second.second);
            }
        }
    }
}

size_t ConnectGraphBuilder::build(Graph* graph) const {
    size_t pair_kmer_num = 0;

    LOG4CXX_DEBUG(logger, boost::format("build graph begin"));
    BOOST_FOREACH(const PairRead& pair_read, _pair_reads) {
        if (pair_read.set1.size() < _K || pair_read.set2.size() < _K) {
            continue;
        }
        pair_kmer_num += addEdge(pair_read.set1, make_complement_dna(pair_read.set2), graph);
        pair_kmer_num += addEdge(pair_read.set2, make_complement_dna(pair_read.set1), graph);
    }
    LOG4CXX_DEBUG(logger, boost::format("build graph end"));

    return pair_kmer_num;
}

struct PairKmerCmp {
    bool operator()(const std::pair<size_t, size_t>& l, const std::pair<size_t, size_t>& r) {
        if (l.first != r.first) {
            return l.first < r.first;
        }
        return l.second < r.second;
    }
};

size_t ConnectGraphBuilder::addEdge(const std::string& read1, const std::string& read2, Graph* graph) const {
    size_t num = 0;

    std::set< std::pair< size_t, size_t >, PairKmerCmp > find_component;
    for (size_t i = 0, j = _K; j <= read1.size() && j <= read2.size(); ++i,++j) {
        Kmer kmer1 = read1.substr(i, _K), kmer2 = read2.substr(i, _K);

        KmerList::const_iterator left = _hash_tbl.find(kmer1), right = _hash_tbl.find(kmer2);
        
        if (left != _hash_tbl.end() && right != _hash_tbl.end() && left->second.first != right->second.first) {
            ++num;

            LOG4CXX_TRACE(logger, boost::format("pair kmre1=%s %d %d") % kmer1 % left->second.first % left->second.second);
            LOG4CXX_TRACE(logger, boost::format("pair kmre2=%s %d %d") % kmer2 % right->second.first % right->second.second);
        
            long distance = _insert_size + left->second.second - right->second.second;
            long left_len = _components[left->second.first].length();
            long overlap = left_len - distance - _K + 1;
            if (overlap > 0) {
                if (overlap > _insert_size)
                    continue;
                distance = left_len - _K + 1; //can change
            }
            if (find_component.find(std::make_pair(left->second.first, right->second.first)) != find_component.end()) {
                graph->addEdge(left->second.first, right->second.first, distance, 1, 0);
            } else {
                find_component.insert( std::make_pair(left->second.first, right->second.first));
                graph->addEdge(left->second.first, right->second.first, distance, 1, 1);
            }
        }
    }

    return num;
}

std::ostream& operator<<(std::ostream& os, const PairReadList& pair_reads) {
    BOOST_FOREACH(const PairRead& pair_read, pair_reads) {
        os << boost::format("%s\t%s") % pair_read.set1 % pair_read.set2 << std::endl; 
    }
    return os;
}

