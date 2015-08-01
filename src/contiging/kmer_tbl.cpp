#include "kmer_tbl.h"
#include "debruijn_graph.h"
#include "kseq.h"

#include <numeric>

#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.kmer_tbl"));

KmerTable::KmerTable(size_t K, bool do_filter) : _K(K), _do_filter(do_filter) {
    BOOST_ASSERT(_K > 0);
}

KmerTable::~KmerTable() {
}

bool KmerTable::read(std::istream& stream) {
    if (!stream) {
        LOG4CXX_WARN(logger, boost::format("invalid input stream"));
        return false;
    }

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table begin"));

    DNASeqReader reader(stream);

    size_t i = 0;
    DNASeq read;
    while (reader.read(read)) {
        if (read.seq.length() > _K) {
            Kmer kmer(read.seq, 0, _K);
            _hash_tbl[kmer]++;
            for (size_t j = _K; j < read.seq.length(); ++j) {
                kmer.pop();
                kmer.push(read.seq[j]);
                _hash_tbl[kmer]++;
            }
        }
        if (i % 100000 == 0) {
            LOG4CXX_DEBUG(logger, boost::format("#%d reads have bean loaded") % (i));
        }
        ++i;
    }

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table end"));

    return true;
}

void KmerTable::buildDeBruijn(DeBruijnGraph* graph) const {
    LOG4CXX_DEBUG(logger, boost::format("construct de bruijn graph begin"));

    BOOST_ASSERT(graph != NULL);

    for (KmerList::const_iterator it = _hash_tbl.begin(); it != _hash_tbl.end(); ++it) {
        graph->addKmer(it->first, it->second);
    }

    LOG4CXX_DEBUG(logger, boost::format("construct de bruijn graph end"));
}

struct Statistics {
    std::pair< double, double > operator()(const std::pair< double, double >& l, const std::pair< Kmer, size_t >& r) const {
        return std::make_pair(l.first + r.second, l.second + r.second * r.second);
    }
};

struct CoverageSquare {
    double operator()(double l, size_t r) const {
        return l + r * r;
    }
};

void KmerTable::statistics(double* average, double* variance) const {
    if (_hash_tbl.size() > 0) {
       	std::pair< double, double > sigma = std::accumulate(_hash_tbl.begin(), _hash_tbl.end(), std::make_pair(0.0, 0.0), Statistics());
        if (average != NULL) {
            *average = sigma.first / _hash_tbl.size();
        }
        if (variance != NULL) {
        }
    }
}
