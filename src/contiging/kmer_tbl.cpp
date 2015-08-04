#include "kmer_tbl.h"
#include "debruijn_graph.h"
#include "kseq.h"

#include <fstream>
#include <numeric>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.kmer_tbl"));

KmerTable::KmerTable(size_t K, double avg_quality, double min_quality, double percent, size_t read_cutoff, bool do_reversed) : _K(K), _avg_quality(avg_quality), _min_quality(min_quality), _percent(percent), _read_cutoff(read_cutoff), _do_reversed(do_reversed) {
    BOOST_ASSERT(_K > 0);
    BOOST_ASSERT(_K <= _read_cutoff);
    BOOST_ASSERT(0 < _percent && _percent <= 1.0);
}

KmerTable::~KmerTable() {
}

bool KmerTable::read(std::istream& stream) {
    if (!stream) {
        LOG4CXX_WARN(logger, boost::format("invalid input stream"));
        return false;
    }

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table begin"));

    DNASeqReader reader(stream, _percent, _read_cutoff);
    DNASeq read;

    while (reader.read(read)) {
        LOG4CXX_TRACE(logger, boost::format("read: %s") % read.seq);

        if (read.seq.length() >= _K) {
            addRead(read);
            if (_do_reversed) {
                // pair-end
                read.make_complement();
                addRead(read);
            }
        }
    }

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table end"));

    return true;
}

bool KmerTable::read(const std::string& file) {
    std::ifstream stream(file.c_str());
    return read(stream);
}

bool KmerTable::read(const std::vector< std::string >& filelist) {
    BOOST_FOREACH(const std::string& file, filelist) {
        if (!read(file)) {
            return false;
        }
    }
    return true;
}

void KmerTable::addRead(const DNASeq& read) {
    Kmer kmer(read.seq, 0, _K);
    if (isValid(read, 0, _K)) {
        _hash_tbl[kmer]++;
    }

    LOG4CXX_TRACE(logger, boost::format("kmer: %s") % kmer);

    for (size_t i = 1,j = _K; j < read.seq.length(); ++i,++j) {
        kmer.pop();
        kmer.push(read.seq[j]);
        if (isValid(read, i, j + 1)) {
            _hash_tbl[kmer]++;
        }

        LOG4CXX_TRACE(logger, boost::format("kmer: %s") % kmer);
    }
}

bool KmerTable::isValid(const DNASeq& read, size_t i, size_t j) const {
    if (_avg_quality > 0 || _min_quality > 0) {
        size_t avg_quality = 0, min_quality = -1;
        for (size_t k = i; k < j; ++k) {
            avg_quality += read.quality[k];
            min_quality = std::min(min_quality, (size_t)read.quality[k]);
        }
        if (j > i) {
            avg_quality /= (j - i);
        }
        if (avg_quality < _avg_quality || min_quality < _min_quality) {
            return false;
        }
    }
    return true;
}

void KmerTable::buildDeBruijn(DeBruijnGraph* graph) const {
    LOG4CXX_DEBUG(logger, boost::format("construct de bruijn graph begin"));

    BOOST_ASSERT(graph != NULL);

    for (KmerList::const_iterator it = _hash_tbl.begin(); it != _hash_tbl.end(); ++it) {
        if (it->second > 1) {
            graph->addKmer(it->first, it->second);
        }
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
        LOG4CXX_DEBUG(logger, boost::format("statistics: count=[%d], sigma=[%f], delta=[%f]") % _hash_tbl.size() % sigma.first % sigma.second);
        if (average != NULL) {
            *average = sigma.first / _hash_tbl.size();
        }
        if (variance != NULL) {
        }
    }
}
