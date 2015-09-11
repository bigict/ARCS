#include "scaffolding.h"
#include "component.h"
#include "contigs.h"
#include "kmer.h"

#include <iostream>

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Scaffolding"));
Scaffolding Scaffolding::_runner;

template< size_t K >
int _Scaffolding_run_(size_t L, const std::string& rootdir, const Properties& options) {
    // check root dir
    boost::filesystem::path root(rootdir);
    if (!boost::filesystem::exists(rootdir) && !boost::filesystem::create_directory(rootdir)) {
        LOG4CXX_ERROR(logger, boost::format("failed to create directory: %s") % rootdir);
        return 1;
    }

    Kmer< K >::length(L); // IMPORTANT: set kmer length

    // read contig file
    ContigList contigs;
    if (!ReadContigs(options.get< std::string >("C"), contigs)) {
        LOG4CXX_ERROR(logger, boost::format("faild to read contigs from file: %s") % options.get< std::string >("C"));
        return 1;
    }
    size_t GENOME_LEN = GenomeLength(contigs, K);
    LOG4CXX_INFO(logger, boost::format("GENOME_LEN=%d") % GENOME_LEN);

    // read component file
    ComponentList components;
    if (!ReadComponents(options.get< std::string >("f"), contigs, K, components)) {
        LOG4CXX_ERROR(logger, boost::format("faild to read components from file: %s") % options.get< std::string >("f"));
        return 1;
    }
    LOG4CXX_INFO(logger, boost::format("Components size=%d") % components.size());

    return 0;
}

int Scaffolding::run(const Properties& options) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    LOG4CXX_DEBUG(logger, "scaffolding begin");

    size_t K = options.get< size_t >("K", 31);
    std::string rootdir = options.get< std::string >("d", ".");

    // process
    if (0 < K && K <= 32) {
        r = _Scaffolding_run_< 32 >(K, rootdir, options);
    } else if (32 < K && K <= 64) {
        r = _Scaffolding_run_< 64 >(K, rootdir, options);
    } else if (64 < K && K <= 96) {
        r = _Scaffolding_run_< 64 >(K, rootdir, options);
    } else if (96 < K && K <= 128) {
        r = _Scaffolding_run_< 128 >(K, rootdir, options);
    } else if (96 < K && K <= 160) {
        r = _Scaffolding_run_< 160 >(K, rootdir, options);
    } else if (160 < K && K <= 192) {
        r = _Scaffolding_run_< 192 >(K, rootdir, options);
    }

    LOG4CXX_DEBUG(logger, "scaffolding end");
    return r;
}

Scaffolding::Scaffolding() : Runner("d:K:C:f:e:1:2:L:P:p:i:r:R:c:h") {
    RUNNER_INSTALL("scaffolding", this, "scaffold");
}

int Scaffolding::printHelps() const {
    std::cout << "arcs scaffolding -K [kmer] -i [input] -d [workdir]" << std::endl;
    return 256;
}

int Scaffolding::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
