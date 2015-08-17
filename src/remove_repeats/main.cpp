#include "uniq_edge_graph.h"

#include <iostream>
#include <string>
#include <vector>

#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

typedef boost::property_tree::ptree Properties;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("remove_repeats.main"));

class RepeatRemover {
public:
    int run(const Properties& options) {
        if (checkOptions(options) != 0) {
            return 1;
        }

	    size_t EDGE_NUM = options.get< size_t >("EDGE_CLUSTER_NUM");

        LOG4CXX_INFO(logger, boost::format("remove repeats begin"));
        LOG4CXX_INFO(logger, boost::format("Edge num=%d") % EDGE_NUM);

        UniqEdgeGraph graph(options.get< size_t >("K"), EDGE_NUM, options.get< size_t >("MAX_OVERLAP"), options.get< size_t >("ITERATION"));
        
        if (!graph.input_edge_length() || !graph.input_edge_position() || !graph.input_edge_link("contig_arc_graph_after_remove_ambigous_arcs") || !graph.input_component()) {
            LOG4CXX_ERROR(logger, boost::format("load data failed."));
            return 2;
        }
    
        graph.linearize();
        
        LOG4CXX_INFO(logger, boost::format("remove repeats end"));

        return 0;
    }

private:
    int checkOptions(const Properties& options);
    int printHelps() const;
};

int RepeatRemover::checkOptions(const Properties& options) {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }

    std::vector< std::string > necessary = boost::assign::list_of("K")("MAX_OVERLAP");
    BOOST_FOREACH(const std::string& c, necessary) {
        if (options.find(c) == options.not_found()) {
            std::cerr << boost::format("The argument -%s is necessary! type -h for more help") % c << std::endl;
            return 1;
        }
    }
    
    std::vector< std::string > intopt = boost::assign::list_of("K");
    BOOST_FOREACH(const std::string& c, intopt) {
        try {
            size_t K = options.get< size_t >(c);
        } catch (std::exception& e) {
            std::cerr << boost::format("The argument -%s should be a integer. type -h for more help") % c << std::endl;
            return 1;
        }
    }

    return 0;
}

int RepeatRemover::printHelps() const {
     std::cerr << "ARCS - version v0.9" << std::endl
         << "Usage:" << std::endl
         << "Arguments:" << std::endl
         << "-K:     select your kmer size." << std::endl
         << "-O:     select your maximum overlap." << std::endl
         << "-i:     which iteration." << std::endl
         << "-s:     select your configure file name." << std::endl
         << "-c:     select a log config file." << std::endl
         << "-h:     help" << std::endl;

    return 1;
}

int main(int argc, char **argv) {
    Properties options;
    {
        // command line options
        Properties cmd;
        std::map< std::string, std::string > alphabet = boost::assign::map_list_of
                ("K", "K")
                ("O", "MAX_OVERLAP")
                ("i", "ITERATION")
                ;
        // command line options
        std::string opt_string("K:d:O:i:s:c:h");
        int opt = -1;
        while ((opt = getopt(argc, argv, opt_string.c_str())) != -1) {
            std::string key(1, (char)opt);
            if (alphabet.find(key) != alphabet.end()) {
                key = alphabet[key];
            }
            if (optarg == NULL) {
                cmd.put(key, NULL);
            } else {
                cmd.put(key, optarg);
            }
        }

        // config log4cxx.
        if (cmd.find("c") != cmd.not_found()) {
            log4cxx::PropertyConfigurator::configure(options.get< std::string >("c"));
        } else {
            log4cxx::BasicConfigurator::configure();
        }

        // load ini options
        if (cmd.find("s") != cmd.not_found()) {
            try {
                boost::property_tree::read_ini(cmd.get< std::string >("s"), options);
            } catch (const boost::property_tree::ini_parser_error& e) {
                LOG4CXX_ERROR(logger, boost::format("load test.cfg failed(%s).") % e.what());
                return 1;
            }
        }

        // merge options
        for (Properties::const_iterator it = cmd.begin(); it != cmd.end(); it++){
            options.put(it->first,it->second.data());
        }
    }
    
    RepeatRemover rr;
    return rr.run(options);
}
