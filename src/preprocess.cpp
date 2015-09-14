#include "preprocess.h"
#include "kmer_dataset.h"

#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Preprocess"));
Preprocess Preprocess::_runner;

class ReadQuality {
public:
    ReadQuality(size_t max_count = -1) : _max_count(max_count) {
    }

    size_t quality(const std::vector< std::string >& filelist, size_t* count) const {
        size_t q = 0, n = 0;

        BOOST_FOREACH(const std::string& file, filelist) {
            std::ifstream stream (file.c_str());
            q += quality(stream, &n);
        }

        if (count != NULL) {
            *count += n;
        }
        return q;
    }
    size_t quality(std::istream& stream, size_t* count) const {
        size_t q = 0, n = 0;

        if (stream) {
            DNASeqReader reader(stream);
            DNASeq read;
            
            while (reader.read(read) && n < _max_count) {
                q += std::accumulate(read.quality.begin(), read.quality.end(), 0);
                n += read.seq.length();
            }
        }

        if (count != NULL) {
            *count += n;
        }
        return q;
    }

    double threshold(const std::vector< std::string >& filelist, double* min_threshold) const {
        size_t n = 0;
        size_t q = quality(filelist, &n);
        return threshold(q, n, min_threshold);
    }

    double threshold(std::istream& stream, double* min_threshold) const {
        size_t n = 0;
        size_t q = quality(stream, &n);
        return threshold(q, n, min_threshold);
    }
    
    double threshold(size_t quality, size_t count, double* min_threshold) const {
        if (count > 0) {
            if (min_threshold != NULL) {
                *min_threshold = (quality / count) * 0.8;
            }
            return (quality / count) * 0.9;
        }
        return 0;
    }

private:
    size_t _max_count;
};

template< size_t K >
int _Preprocess_run_(size_t L, size_t buckets, double avg_quality, double min_quality, double percent, size_t read_cutoff, bool do_reverse, std::istream& is, std::ostream& os) {
    Kmer< K >::length(L); // IMPORTANT: set kmer length
    KmerDataSet< K > tbl(buckets, avg_quality, min_quality, percent, read_cutoff, do_reverse);
    if (!tbl.read(is)) {
        return 1;
    }
    if (!tbl.write(os)) {
        return 1;
    }
    return 0;
}

int Preprocess::run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    LOG4CXX_DEBUG(logger, "preprocess reads begin");

    // parameters
    size_t K = options.get< size_t >("K", 31);
    size_t buckets = options.get< size_t >("n", 0);
    double avg_quality = 0, min_quality = 0, percent = 1.0;
    std::vector< std::string > filelist;
    {
        std::string input = options.get< std::string >("i", "");
        boost::algorithm::split(filelist, input, boost::algorithm::is_any_of(":"));
        std::cout << boost::algorithm::join(filelist, "$")<< std::endl;
    }
    if (options.find("E") != options.not_found() && !filelist.empty()) {
        ReadQuality statistics;
        avg_quality = statistics.threshold(filelist, &min_quality);
        percent = 0.95;
    }
    size_t read_cutoff = options.get< size_t >("READ_LENGTH_CUTOFF", -1);
    bool do_reverse = options.find("S") == options.not_found();

    // input & output opened
    std::istream* is = &std::cin;
    if (false) {
        std::string file = options.get< std::string >("i");
        is = new std::ifstream(file.c_str());

        LOG4CXX_DEBUG(logger, boost::format("input: %s") % file);
    }
    if (is == &std::cin) {
        std::cin.sync_with_stdio(false);
    }
    std::ostream* os = &std::cout;
    if (options.find("o") != options.not_found()) {
        std::string file = options.get< std::string >("o");
        os = new std::ofstream(file.c_str());

        LOG4CXX_DEBUG(logger, boost::format("output: %s") % file);
    }

    LOG4CXX_DEBUG(logger, boost::format("parameters: K=[%d], buckets=[%d], avg_quality=[%lf],min_quality=[%lf],read_cutoff=[%d],do_reverse=[%d]") % K % buckets % avg_quality % min_quality % read_cutoff % do_reverse);

    // process
    if (0 < K && K <= 32) {
        r = _Preprocess_run_< 32 >(K, buckets, avg_quality, min_quality, percent, read_cutoff, do_reverse, *is, *os);
    } else if (32 < K && K <= 64) {
        r = _Preprocess_run_< 64 >(K, buckets, avg_quality, min_quality, percent, read_cutoff, do_reverse, *is, *os);
    } else if (64 < K && K <= 96) {
        r = _Preprocess_run_< 64 >(K, buckets, avg_quality, min_quality, percent, read_cutoff, do_reverse, *is, *os);
    } else if (96 < K && K <= 128) {
        r = _Preprocess_run_< 128 >(K, buckets, avg_quality, min_quality, percent, read_cutoff, do_reverse, *is, *os);
    } else if (96 < K && K <= 160) {
        r = _Preprocess_run_< 160 >(K, buckets, avg_quality, min_quality, percent, read_cutoff, do_reverse, *is, *os);
    } else if (160 < K && K <= 192) {
        r = _Preprocess_run_< 192 >(K, buckets, avg_quality, min_quality, percent, read_cutoff, do_reverse, *is, *os);
    } else {
        r = 1;
    }

    // input & output closed
    if (os != &std::cout) {
        delete os;
    }
    if (is != &std::cin) {
        delete is;
    }

    LOG4CXX_DEBUG(logger, "preprocess reads end");
    return r;
}

Preprocess::Preprocess() : Runner("c:s:K:n:i:o:E:h", boost::assign::map_list_of('E', "READ_LENGTH_CUTOFF")) {
    RUNNER_INSTALL("preprocess", this, "preprocess");
}

int Preprocess::printHelps() const {
    std::cout << "arcs preprocess -K [kmer] -i [input] -o [output] -E [read_cutoff] -n [buckets]" << std::endl;
    return 256;
}

int Preprocess::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
