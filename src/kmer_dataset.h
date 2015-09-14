#ifndef kmer_dataset_h_
#define kmer_dataset_h_

#include "kmer.h"
#include "kseq.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/assert.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

template< size_t K >
class KmerDataSet {
public:
    KmerDataSet(size_t n, double avg_quality=0, double min_quality=0, double percent=1.0, size_t read_cutoff=-1, bool do_reversed=true) : _hash_tbl(n), _avg_quality(avg_quality), _min_quality(min_quality), _percent(percent), _read_cutoff(read_cutoff), _do_reversed(do_reversed) {
        BOOST_ASSERT(K > 0);
        BOOST_ASSERT(K <= _read_cutoff);
        BOOST_ASSERT(0 < _percent && _percent <= 1.0);
        
        LOG4CXX_DEBUG(logger, boost::format("buckets=[%d],avg_quality=[%lf],min_quality=[%lf],percent=[%lf],read_cutoff=[%d],do_reversed=[%d]") % n % avg_quality % min_quality % percent % read_cutoff % do_reversed);
    }

    bool read(std::istream& stream) {
        if (!stream) {
            LOG4CXX_WARN(logger, boost::format("invalid input stream"));
            return false;
        }

        LOG4CXX_DEBUG(logger, boost::format("construct kmer dataset from stream begin"));

        DNASeqReader reader(stream, _percent, _read_cutoff);
        DNASeq read;

        while (reader.read(read)) {
            LOG4CXX_TRACE(logger, boost::format("read: %s") % read.seq);

            if (read.seq.length() >= Kmer< K >::length()) {
                addRead(read);
                if (_do_reversed) {
                    // pair-end
                    read.make_complement();
                    addRead(read);
                }
            }
        }

        LOG4CXX_DEBUG(logger, boost::format("construct kmer dataset from stream end"));

        return true;
    }
    bool read(const std::string& file) {
        std::ifstream stream(file.c_str());
        return read(stream);
    }
    bool read(const std::vector< std::string >& filelist) {
        BOOST_FOREACH(const std::string& file, filelist) {
            if (!read(file)) {
                return false;
            }
        }
        return true;
    }
    bool write(std::ostream& stream) const {
        double mean = 0, var = 0;
        statistics(&mean, &var);

        stream << boost::format("@attribute mean %f") % mean << std::endl;
        stream << boost::format("@attribute variance %f") % var << std::endl;
        stream << boost::format("@data") << std::endl;
        for (typename KmerList::const_iterator it = _hash_tbl.begin(); it != _hash_tbl.end(); ++it) {
            if (it->second > 1) {
                stream << boost::format("%s\t%d") % it->first % it->second << std::endl;
            }
        }
        return true;
    }

    size_t size() const {
        return _hash_tbl.size();
    }

    void statistics(double* average, double* variance) const;

private:
    void addRead(const DNASeq& read) {
        size_t L = Kmer< K >::length();

        Kmer< K > kmer(read.seq, 0, L);
        if (isValid(read, 0, L)) {
            //_hash_tbl[kmer]++;
            typename KmerList::iterator k = _hash_tbl.find(kmer);
            if (k != _hash_tbl.end()) {
                ++k->second;
            } else {
                _hash_tbl.insert(std::make_pair(kmer, 1));
            }
        }

        LOG4CXX_TRACE(logger, boost::format("kmer: %s") % kmer);

        for (size_t i = 1,j = L; j < read.seq.length(); ++i,++j) {
            kmer.shift(read.seq[j]);
            if (isValid(read, i, j + 1)) {
                //_hash_tbl[kmer]++;
                typename KmerList::iterator k = _hash_tbl.find(kmer);
                if (k != _hash_tbl.end()) {
                    ++k->second;
                } else {
                    _hash_tbl.insert(std::make_pair(kmer, 1));
                }
            }

            LOG4CXX_TRACE(logger, boost::format("kmer: %s") % kmer);
        }
    }
    bool isValid(const DNASeq& read, size_t i, size_t j) const {
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

    typedef KmerTable< K, size_t > KmerList; 
    KmerList _hash_tbl; 
    double _avg_quality;
    double _min_quality;
    double _percent;
    size_t _read_cutoff;
    bool _do_reversed;

    static log4cxx::LoggerPtr logger;
};

template< size_t K >
log4cxx::LoggerPtr KmerDataSet< K >::logger(log4cxx::Logger::getLogger("arcs.KmerDataSet"));

typedef boost::accumulators::accumulator_set< size_t, boost::accumulators::stats< boost::accumulators::tag::count, boost::accumulators::tag::mean, boost::accumulators::tag::moment< 2 > > > _Accumulator;

template< size_t K >
struct _Statistics {
    _Statistics(size_t min_threshold, _Accumulator& acc) : _min_threshold(min_threshold), _acc(acc) {
    }
    void operator()(const std::pair< Kmer< K >, size_t >& item) const {
        if (item.second >= _min_threshold) {
            _acc(item.second);
        }
    }
private:
    _Accumulator& _acc;
    size_t _min_threshold;
};

template< size_t K >
void KmerDataSet< K >::statistics(double* average, double* variance) const {
    if (_hash_tbl.size() > 0) {
        _Accumulator acc;
        std::for_each(_hash_tbl.begin(), _hash_tbl.end(), _Statistics< K >(4, acc));

        LOG4CXX_DEBUG(logger, boost::format("statistics: tbl_size=[%d], mean=[%f], moment<2>=[%f],count=[%d]") % _hash_tbl.size() % boost::accumulators::mean(acc) % boost::accumulators::moment< 2 >(acc) % boost::accumulators::count(acc));

        if (average != NULL) {
            *average = boost::accumulators::mean(acc);
        }
        if (variance != NULL) {
            *variance = std::sqrt(boost::accumulators::moment< 2 >(acc) - std::pow(boost::accumulators::mean(acc), 2));
        }
    }
}

#endif // kmer_dataset_h_
