#include "kmer_tbl.h"
#include "debruijn_graph.h"
#include "kseq.h"

#include <fstream>
#include <numeric>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/bind.hpp>
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

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table from stream begin"));

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

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table from stream end"));

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

bool KmerTable::write(std::ostream& stream) const {
    for (KmerList::const_iterator it = _hash_tbl.begin(); it != _hash_tbl.end(); ++it) {
        if (it->second > 1) {
            stream << boost::format("%s\t%d") % it->first % it->second << std::endl;
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
    LOG4CXX_DEBUG(logger, boost::format("construct de bruijn graph from kmer_tbl begin"));

    BOOST_ASSERT(graph != NULL);

    for (KmerList::const_iterator it = _hash_tbl.begin(); it != _hash_tbl.end(); ++it) {
        if (it->second > 1) {
            graph->addKmer(it->first, it->second);
        }
    }

    LOG4CXX_DEBUG(logger, boost::format("construct de bruijn graph kmer_tbl end"));
}

typedef boost::accumulators::accumulator_set< size_t, boost::accumulators::stats< boost::accumulators::tag::count, boost::accumulators::tag::mean, boost::accumulators::tag::moment< 2 > > > Accumulator;
struct Statistics {
    Statistics(size_t min_threshold, Accumulator& acc) : _min_threshold(min_threshold), _acc(acc) {
    }
    void operator()(const std::pair< Kmer, size_t >& item) const {
        if (item.second >= _min_threshold) {
            _acc(item.second);
        }
    }
private:
    Accumulator& _acc;
	size_t _min_threshold;
};

void KmerTable::statistics(double* average, double* variance) const {
    if (_hash_tbl.size() > 0) {
        Accumulator acc;
        std::for_each(_hash_tbl.begin(), _hash_tbl.end(), Statistics(4, acc));

        LOG4CXX_DEBUG(logger, boost::format("statistics: tbl_size=[%d], mean=[%f], moment<2>=[%f],count=[%d]") % _hash_tbl.size() % boost::accumulators::mean(acc) % boost::accumulators::moment< 2 >(acc) % boost::accumulators::count(acc));

        if (average != NULL) {
            *average = boost::accumulators::mean(acc);
        }
        if (variance != NULL) {
            *variance = std::sqrt(boost::accumulators::moment< 2 >(acc) - std::pow(boost::accumulators::mean(acc), 2));
        }
    }
}
