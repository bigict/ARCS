#include "contiging.h"
#include "debruijn_graph.h"

#include <iostream>

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Contiging"));
Contiging Contiging::_runner;

template< size_t K >
int _Contiging_run_(size_t L, size_t loops, std::istream& is, const std::string& rootdir) {
    // check root dir
    boost::filesystem::path root(rootdir);
    if (!boost::filesystem::exists(rootdir) && !boost::filesystem::create_directory(rootdir)) {
        LOG4CXX_ERROR(logger, boost::format("failed to create directory: %s") % rootdir);
        return 1;
    }

    Kmer< K >::length(L); // IMPORTANT: set kmer length
    DeBruijnGraph< K > g;
    
    // load
    if (!g.read(is)) {
        return 1;
    }

    // output
    {
        boost::filesystem::ofstream stream(rootdir / boost::filesystem::path("condensed_de_bruijn_graph_before_trimming.data"));
        if (!stream) {
            return 1;
        }
        stream << g;
    }

    // trimming
    for (size_t i = 0; i < loops; ++i) {
        g.removeNoise();
    }

    // output
    {
        boost::filesystem::ofstream stream(rootdir / boost::filesystem::path("condensed_de_bruijn_graph_after_trimming.data"));
        if (!stream) {
            return 1;
        }
        stream << g;
    }
    // parameters
    {
        boost::filesystem::ofstream stream(rootdir / boost::filesystem::path("contig_parameter"));
        if (!stream) {
            return 1;
        }
        stream << boost::format("lambda=%f") % g.average() << std::endl;
    }

    return 0;
}

int Contiging::run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    LOG4CXX_DEBUG(logger, "contiging begin");

    size_t K = options.get< size_t >("K", 31);
    size_t loops = options.get< size_t >("l", 2);
    std::string rootdir = options.get< std::string >("d", ".");

    // input & output opened
    std::istream* is = &std::cin;
    if (options.find("i") != options.not_found()) {
        std::string file = options.get< std::string >("i");
        is = new std::ifstream(file.c_str());

        LOG4CXX_DEBUG(logger, boost::format("input: %s") % file);
    }
    if (is == &std::cin) {
        std::cin.sync_with_stdio(false);
    }

    LOG4CXX_DEBUG(logger, boost::format("parameters: K=[%d],workdir=[%s]") % K % rootdir);

    K = K - 1; // (K-1)mer
    // process
    if (0 < K && K <= 32) {
        r = _Contiging_run_< 32 >(K, loops, *is, rootdir);
    } else if (32 < K && K <= 64) {
        r = _Contiging_run_< 64 >(K, loops, *is, rootdir);
    } else if (64 < K && K <= 96) {
        r = _Contiging_run_< 64 >(K, loops, *is, rootdir);
    } else if (96 < K && K <= 128) {
        r = _Contiging_run_< 128 >(K, loops, *is, rootdir);
    } else if (96 < K && K <= 160) {
        r = _Contiging_run_< 160 >(K, loops, *is, rootdir);
    } else if (160 < K && K <= 192) {
        r = _Contiging_run_< 192 >(K, loops, *is, rootdir);
    }

    // input & output closed
    if (is != &std::cin) {
        delete is;
    }

    LOG4CXX_DEBUG(logger, "contiging end");
    return r;
}

Contiging::Contiging() : Runner("c:s:K:i:d:h") {
    RUNNER_INSTALL("contiging", this, "contigs");
}

int Contiging::printHelps() const {
    std::cout << "arcs contiging -K [kmer] -i [input] -d [workdir]" << std::endl;
    return 256;
}

int Contiging::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
