#include "kmer_tbl.h"
#include "debruijn_graph.h"
#include "kseq.h"

#include <fstream>
#include <numeric>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.kmer_tbl"));

KmerTable::KmerTable(size_t K, bool do_filter, bool do_reversed) : _K(K), _do_filter(do_filter), _do_reversed(do_reversed) {
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
    DNASeq read;

    double average_threshold = 0, min_threshold = 0;

    if (_do_filter) {
        size_t pos = stream.tellg();
        size_t count = 0, quality = 0;

        while (reader.read(read)) {
            quality += std::accumulate(read.quality.begin(), read.quality.end(), 0);
            count += read.seq.length();
        }
        stream.clear();
        stream.seekg(pos, stream.beg);

        if (count > 0) {
            average_threshold =  (quality / count) * 0.9;
            min_threshold = (quality / count) * 0.8;
        }
    }

    while (reader.read(read)) {
        LOG4CXX_TRACE(logger, boost::format("read: %s") % read.seq);

        if (read.seq.length() > _K) {
            addRead(read);
            if (_do_reversed) {
                // pair-end
                std::reverse(read.seq.begin(), read.seq.end());
                std::reverse(read.quality.begin(), read.quality.end());
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
    _hash_tbl[kmer]++;

    LOG4CXX_TRACE(logger, boost::format("kmer: %s") % kmer);

    for (size_t j = _K; j < read.seq.length(); ++j) {
        kmer.pop();
        kmer.push(read.seq[j]);
        _hash_tbl[kmer]++;

        LOG4CXX_TRACE(logger, boost::format("kmer: %s") % kmer);
    }
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
        LOG4CXX_DEBUG(logger, boost::format("statistics: count=[%d], sigma=[%f], delta=[%f]") % _hash_tbl.size() % sigma.first % sigma.second);
        if (average != NULL) {
            *average = sigma.first / _hash_tbl.size();
        }
        if (variance != NULL) {
        }
    }
}
