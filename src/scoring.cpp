#include "scoring.h"
#include "uniq_edge_graph.h"

#include <iostream>
#include <tr1/tuple>

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Scoring"));
Scoring Scoring::_runner;

int Scoring::run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }
    boost::filesystem::path workdir(options.get< std::string >("d", "."));
    if (!boost::filesystem::exists(workdir) && !boost::filesystem::create_directory(workdir)) {
        LOG4CXX_ERROR(logger, boost::format("failed to create directory: %s") % workdir);
        return 1;
    }


    LOG4CXX_DEBUG(logger, "scoring begin");

    size_t K = options.get< size_t >("K");
    size_t ITERATION = options.get< size_t >("ITERATION");
    size_t MAX_OVERLAP = options.get< size_t >("MAX_OVERLAP");

    LOG4CXX_DEBUG(logger, boost::format("parameters: K=[%d],ITERATION=[%d],MAX_OVERLAP=[%d]") % K % ITERATION % MAX_OVERLAP);

    UniqEdgeGraph graph(K, MAX_OVERLAP, ITERATION);
    // load data
    {
        typedef bool(UniqEdgeGraph::*LoadDataPtr)(std::istream&);
        typedef std::tr1::tuple< std::string, LoadDataPtr> File2FuncPtr;
        std::vector< File2FuncPtr > file2func_list = boost::assign::list_of
            (std::tr1::make_tuple("edge_cluster_len_%ld", &UniqEdgeGraph::input_edge_length))
            (std::tr1::make_tuple("edge_cluster_pos_%ld", &UniqEdgeGraph::input_edge_position))
            (std::tr1::make_tuple("contig_arc_graph_after_remove_ambigous_arcs_%ld", &UniqEdgeGraph::input_edge_link))
            (std::tr1::make_tuple("component_%ld", &UniqEdgeGraph::input_component))
            ;

        BOOST_FOREACH(const File2FuncPtr& file2func, file2func_list) {
            std::string file = boost::str(boost::format(std::tr1::get< 0 >(file2func)) % ITERATION);
            boost::filesystem::ifstream stream(workdir / file);
            if (!boost::bind(std::tr1::get< 1 >(file2func), &graph, _1)(stream)) {
                LOG4CXX_ERROR(logger, boost::format("failed to load %s") % file);
                return 2;
            }
        }
    }
    // resolve conflicts && output
    {
        boost::filesystem::ofstream stream(workdir / boost::filesystem::path(
                    boost::str(boost::format("component_%ld") % (ITERATION + 1))
            ));
        if (!stream) {
            return 1;
        }
        graph.linearize(stream);
    }

    LOG4CXX_DEBUG(logger, "scoring end");
    return 0;
}

Scoring::Scoring() : Runner("c:s:K:d:O:i:h", boost::assign::map_list_of('O', "MAX_OVERLAP")('i', "ITERATION")) {
    RUNNER_INSTALL("scoring", this, "scoring");
}

int Scoring::printHelps() const {
    std::cout << "arcs scoring -K [kmer] -O [max overlap] -i [iteration] -d [workdir]" << std::endl;
    return 256;
}

int Scoring::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
