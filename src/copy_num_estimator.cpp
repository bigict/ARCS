#include "copy_num_estimator.h"
#include "constant.h"
#include "condensed_debruijn_graph_reader.h"
#include "contigs.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.CopyNumEstimator"));

int CopyNumEstimator::run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    double lambda = options.get< double >("lambda");

    LOG4CXX_DEBUG(logger, "copy_num_estimator begin");

    ContigList contigs;
    // read contigs
    {
        std::istream* is = &std::cin;
        if (options.find("i") != options.not_found()) {
            std::string file = options.get< std::string >("i");
            is = new std::ifstream(file.c_str());
        }
        if (is == &std::cin) {
            std::cin.sync_with_stdio(false);
        }

        size_t edge_num = 0;
        CondensedDeBruijnGraphReader reader(*is);
        CondensedDeBruijnGraphEdge edge;
        while (reader.read(edge)) {
            int copy_num = (int)std::round(edge.coverage / lambda);
            if (copy_num >= 0) {
                contigs.push_back(Contig(edge.seq, copy_num));
            }
            ++edge_num;
        }

        LOG4CXX_DEBUG(logger, boost::format("Number of condensed edges = %d") % edge_num);
        LOG4CXX_DEBUG(logger, boost::format("Number of Contigs = %d") % contigs.size());

        if (is != &std::cin) {
            delete is;
        }
    }
    // write contigs
    {
        std::string file = options.get< std::string >("G");
        if (!WriteContigs(file, contigs)) {
            return 1;
        }
    }
    // write components
    {
        std::string file = options.get< std::string >("C");
        std::ofstream stream(file.c_str());
        if (!stream) {
            return 1;
        }

        size_t contig_idx = 0, component_idx = 0;
        BOOST_FOREACH(const Contig& contig, contigs) {
            if (contig.copy_num <= 1) {
                stream << boost::format(">component\t%d\n%d\n\n") % component_idx % contig_idx;
                ++component_idx;
            }
            ++contig_idx;
        }

        LOG4CXX_DEBUG(logger, boost::format("Number of Uniques = %d") % component_idx);
    }

    LOG4CXX_DEBUG(logger, "copy_num_estimator end");
    return r;
}

CopyNumEstimator CopyNumEstimator::_runner;

CopyNumEstimator::CopyNumEstimator() : Runner("c:s:i:G:C:h") {
    RUNNER_INSTALL("copy_num_estimate", this, "estimate copy number of contigs");
}

int CopyNumEstimator::printHelps() const {
    std::cout << "arcs copy_num_estimate -i [input] -G [cdbg] -C [component]" << std::endl;
    return 256;
}

int CopyNumEstimator::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
