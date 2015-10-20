#ifndef pair_read_h_
#define pair_read_h_

#include "component.h"
#include "kmer_tbl.h"
#include "kseq.h"

#include <iostream>
#include <string>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <log4cxx/logger.h>

struct PairRead {
    std::string set1;
    std::string set2;
};

typedef std::vector< PairRead > PairReadList;
bool ReadPairReads(std::istream& stream1, std::istream& stream2, PairReadList& pair_reads);
bool ReadPairReads(const std::string& file1, const std::string& file2, PairReadList& pair_reads);

template< size_t K >
class InsertSizeEstimater {
public:
    InsertSizeEstimater(size_t L, size_t insert_size, const PairReadList& pair_reads, const KmerTable< K, KmerPosition >& hash_tbl) : _K(L), _insert_size(insert_size), _pair_reads(pair_reads), _hash_tbl(hash_tbl) {
    }

    void estimate(size_t* insert_size, double* delta);
private:
    typedef std::vector< int > InsertSizeDistr;
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


    void estimateOnePR(const std::string& read1, const std::string& read2, InsertSizeDistr& insert_size_distr) {
        LOG4CXX_TRACE(logger, boost::format("InsertSizeEstimater::estimateOnePR %s %s") % read1 % read2);

        for (size_t i = 0, j = _K; j <= read1.size() && j <= read2.size(); ++i,++j) {
            Kmer< K > kmer1 = read1.substr(i, _K), kmer2 = read2.substr(i, _K);
            typename KmerTable< K, KmerPosition >::const_iterator left = _hash_tbl.find(kmer1), right = _hash_tbl.find(kmer2);

            if (left != _hash_tbl.end() && right != _hash_tbl.end() && left->second.first == right->second.first) {
                LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre1=%s %d %d") % kmer1 % left->second.first % left->second.second);
                LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre2=%s %d %d") % kmer2 % right->second.first % right->second.second);

                if (right->second.second - left->second.second > _insert_size / 4 
                        && right->second.second - left->second.second < 2 * _insert_size) { // use abs(x) ?????????
                    insert_size_distr.push_back(right->second.second - left->second.second);
                }
            }
        }
    }

    size_t _K;
    size_t _insert_size;
    const PairReadList& _pair_reads;
    const KmerTable< K, KmerPosition >& _hash_tbl;

    static log4cxx::LoggerPtr logger;
};

template< size_t K>
log4cxx::LoggerPtr InsertSizeEstimater< K >::logger(log4cxx::Logger::getLogger("arcs.PairRead"));

template< size_t K >
void InsertSizeEstimater< K >::estimate(size_t* insert_size, double* delta) {
    LOG4CXX_DEBUG(logger, boost::format("pair_reads=[%d]") % _pair_reads.size());

    if (insert_size != NULL) {
        *insert_size = 0;
    }
    if (delta != NULL) {
        *delta = 0.0;
    }

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
    if (boost::accumulators::count(acc) > 0) {
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
}

#endif // pair_read_h_
