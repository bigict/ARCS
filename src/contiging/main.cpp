#include <iostream>
#include <fstream>
#include <vector>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

#include "kmer_tbl.h"
#include "debruijn_graph.h"

typedef boost::property_tree::ptree Properties;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.main"));

class Contiging {
public:
    int run(const Properties& options) {
        if (checkOptions(options) != 0) {
            return 1;
        }

        KmerTable tbl(options.get< size_t >("K"), options.find("E") != options.not_found(), options.find("S") == options.not_found());
        std::vector< std::string > filelist = boost::assign::list_of(options.get< std::string >("q1"))(options.get< std::string >("q2"));
        if (!tbl.read(filelist)) {
            LOG4CXX_ERROR(logger, boost::format("faild to open file: %s") % options.get< std::string >("q1"));
            return -1;
        }

        DeBruijnGraph graph(tbl);
        graph.compact();
        graph.removeNoise();
        
        return 0;
    }
private:
    int checkOptions(const Properties& options);
    int printHelps() const;
};

int Contiging::checkOptions(const Properties& options) {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }

    std::vector< std::string > necessary = boost::assign::list_of("d")("s")("p")("K");
    BOOST_FOREACH(const std::string& c, necessary) {
        if (options.find(c) == options.not_found()) {
            std::cerr << boost::format("The argument -%s is necessary! type -h for more help") % c << std::endl;
            return 1;
        }
    }
    
    try {
        size_t K = options.get< size_t >("K");
    } catch (std::exception& e) {
        std::cerr << boost::format("The argument -K should be a integer. type -h for more help") << std::endl;
    }
    
    return 0;
}

int Contiging::printHelps() const {
     std::cerr << "ARCS - version v0.9" << std::endl
         << "Usage:" << std::endl
         << "Arguments:" << std::endl
         << "-d:     select your work directory." << std::endl
         << "-s:     select your configure file name." << std::endl
         << "-E:     select whether you want to filter out bad reads." << std::endl
         << "-p:     select the number of the threads you want to use. The argument--p must be a figure more than zero." << std::endl
         << "-K:     select your kmer size." << std::endl
         << "-S:     select whether you want to read pair ends." << std::endl
         << "-c:     select a log config file." << std::endl
         << "-h:     help" << std::endl;

    return 1;
}

int main(int argc, char* argv[]) {
    Properties options;
    {
        // command line options
        Properties cmd;
        std::string opt_string("s:d:K:ESp:c:h");
        int opt = -1;
        while ((opt = getopt(argc, argv, opt_string.c_str())) != -1) {
            std::string key(1, (char)opt);
            if (optarg == NULL) {
                cmd.put(key, NULL);
            } else {
                cmd.put(key, optarg);
            }
        }

        // config log4cxx.
        if (cmd.find("c") != cmd.not_found()) {
            log4cxx::PropertyConfigurator::configure(cmd.get< std::string >("c"));
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

    // build contigs.
    Contiging c;
    return c.run(options);
}
