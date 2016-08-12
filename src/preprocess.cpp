#include "preprocess.h"
#include "constant.h"
#include "kmer_dataset.h"

#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Preprocess"));

class ReadQuality {
public:
    ReadQuality(size_t max_count = -1) : _max_count(max_count) {
    }

    size_t quality(const std::vector< std::string >& filelist, size_t* count) const {
        size_t q = 0, n = 0;

        BOOST_FOREACH(const std::string& file, filelist) {
            if(boost::algorithm::ends_with(file, ".gz")) {
                std::ifstream is(file.c_str(), std::ios_base::in | std::ios_base::binary);
                boost::iostreams::filtering_istream stream;
                stream.push(boost::iostreams::gzip_decompressor());
                stream.push(is);
                q += quality(stream, &n);
            } else {
                std::ifstream stream (file.c_str());
                q += quality(stream, &n);
            }
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
int _Preprocess_run_(size_t L, const Properties& options, const Arguments& arguments) {
    size_t r = 0;

    Kmer< K >::length(L); // IMPORTANT: set kmer length

    // parameters
    size_t buckets = options.get< size_t >("n", 0);
    double avg_quality = 0, min_quality = 0, percent = 1.0;
    std::vector< std::string > filelist;
    for(std::string files : arguments) {
        std::vector< std::string > res;
        boost::split(res, files, boost::is_any_of(":"));
        std::copy(res.begin(), res.end(), std::back_inserter(filelist));
    }
    {
        std::string file = boost::algorithm::join(filelist, ":");
        LOG4CXX_DEBUG(logger, boost::format("input: %s") % file);
    }
    if (options.find("E") != options.not_found() && !filelist.empty()) {
        ReadQuality statistics;
        avg_quality = statistics.threshold(filelist, &min_quality);
        percent = 0.95;
    }
    size_t read_cutoff = options.get< size_t >("READ_LENGTH_CUTOFF", -1);
    bool do_reverse = options.find("S") == options.not_found();
    bool filt_kmer_with_N = options.find("F") == options.not_found();

    LOG4CXX_DEBUG(logger, boost::format("parameters: K=[%d], buckets=[%d], avg_quality=[%lf],min_quality=[%lf],read_cutoff=[%d],do_reverse=[%d],filt_kmer_with_N=[%d]") % L % buckets % avg_quality % min_quality % read_cutoff % do_reverse % filt_kmer_with_N);

    // construct dataset
    KmerDataSet< K > tbl(buckets, avg_quality, min_quality, percent, read_cutoff, do_reverse, filt_kmer_with_N);

    // read
    if (r == 0){
        if (!filelist.empty()) {
            r = tbl.read(filelist) ? 0 : 1;
        } else {
            std::cin.sync_with_stdio(false);
            r = tbl.read(std::cin) ? 0 : 1;
        }
    }
    LOG4CXX_DEBUG(logger, boost::format("table=[%d]") % tbl.size());
    // write
    if (r == 0) {
        size_t threshold = options.get< size_t >("THRESHOLD", 1);
        if (options.find("o") != options.not_found()) {
            std::string file = options.get< std::string >("o");
            std::ofstream stream(file.c_str());
            r = tbl.write(stream, threshold) ? 0 : 1;

            LOG4CXX_DEBUG(logger, boost::format("output: %s") % file);
        } else {
            r = tbl.write(std::cout, threshold) ? 0 : 1;
        }
    }

    return r;
}

int Preprocess::run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    LOG4CXX_DEBUG(logger, "preprocess reads begin");

    // parameters
    size_t K = options.get< size_t >("K", kKmerSize);

    // process
    if (0 < K && K <= 32) {
        r = _Preprocess_run_< 32 >(K, options, arguments);
    } else if ( 32 < K && K <= 64) {
        r = _Preprocess_run_< 64 >(K, options, arguments);
    } else if ( 64 < K && K <= 96) {
        r = _Preprocess_run_< 96 >(K, options, arguments);
    } else if ( 96 < K && K <= 128) {
        r = _Preprocess_run_<128 >(K, options, arguments);
    } else if (128 < K && K <= 160) {
        r = _Preprocess_run_<169 >(K, options, arguments);
    } else if (160 < K && K <= 192) {
        r = _Preprocess_run_<192 >(K, options, arguments);
    } else {
        r = 1;
    }

    LOG4CXX_DEBUG(logger, "preprocess reads end");
    return r;
}

Preprocess Preprocess::_runner;

Preprocess::Preprocess() : Runner("c:s:K:n:d:i:o:t:ESe:Fh", boost::assign::map_list_of('e', "READ_LENGTH_CUTOFF")('t', "THRESHOLD")) {
    RUNNER_INSTALL("preprocess", this, "filter and quality-trim reads");
}

int Preprocess::printHelps() const {
    std::cout << "usage: arcs preprocess [arguments] [<inputs>]" << std::endl;
    std::cout << std::endl;
    std::cout << "\t-c[=<file>]    log config file, default ./log4cxx.properties" << std::endl;
    std::cout << "\t-s[=<file>]    assembly config file" << std::endl;
    std::cout << "\t-K[=<number>]  kmer size, default 31" << std::endl;
    std::cout << "\t-n[=<number>]  initial hash bucket size, default 0" << std::endl;
    //std::cout << "\t-d" << endl;
    //std::cout << "\t-i" << endl;
    std::cout << "\t-o[=<file>]    kmer output file, default stdout" << std::endl;
    std::cout << "\t-t[=<number>]  filter out kmers whose coverage <= number, default 1" << std::endl;
    std::cout << "\t-E             filter out low quality kmers" << std::endl;
    std::cout << "\t-S             do not reverse original reads" << std::endl;
    std::cout << "\t-e[=<number>]  read length cutoff, default -1" << std::endl;
    std::cout << "\t-F             filter out kmers with N, default random" << std::endl;
    std::cout << "\t-h             help" << std::endl;
    std::cout << std::endl;
    return 256;
}

int Preprocess::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
