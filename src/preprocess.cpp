#include "preprocess.h"
#include "kmer_dataset.h"

#include <iostream>

#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Preprocess"));
Preprocess Preprocess::_runner;

template< size_t K >
int _Preprocess_run_(size_t buckets, std::istream& is, std::ostream& os) {
    KmerDataSet< K > tbl(buckets);
    if (!tbl.read(is)) {
        return 1;
    }
    if (!tbl.write(os)) {
        return 1;
    }
    return 0;
}

int Preprocess::run(const Properties& options) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    LOG4CXX_DEBUG(logger, "preprocess reads begin");

    size_t K = options.get< size_t >("K", 31);
    size_t buckets = options.get< size_t >("n", 0);

    // input & output opened
    std::istream* is = &std::cin;
    if (options.find("i") != options.not_found()) {
        std::string file = options.get< std::string >("i");
        is = new std::ifstream(file.c_str());
    }
    if (is == &std::cin) {
        std::cin.sync_with_stdio(false);
    }
    std::ostream* os = &std::cout;
    if (options.find("o") != options.not_found()) {
        std::string file = options.get< std::string >("o");
        os = new std::ofstream(file.c_str());
    }

    LOG4CXX_DEBUG(logger, boost::format("parameters: K=[%d], buckets=[%d]") % K % buckets);

    // process
    if (0 < K && K <= 32) {
        r = _Preprocess_run_< 32 >(buckets, *is, *os);
    } else if (32 < K && K <= 64) {
        r = _Preprocess_run_< 64 >(buckets, *is, *os);
    } else if (64 < K && K <= 96) {
        r = _Preprocess_run_< 64 >(buckets, *is, *os);
    } else if (96 < K && K <= 128) {
        r = _Preprocess_run_< 128 >(buckets, *is, *os);
    } else if (96 < K && K <= 160) {
        r = _Preprocess_run_< 160 >(buckets, *is, *os);
    } else if (160 < K && K <= 192) {
        r = _Preprocess_run_< 192 >(buckets, *is, *os);
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

int Preprocess::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}

int Preprocess::printHelps() const {
    std::cout << "arcs preprocess -K [kmer] -i [input] -o [output] -n [buckets]" << std::endl;
    return 256;
}
