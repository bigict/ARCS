#include <iostream>
#include <fstream>

#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

#include "debruijn_graph.h"

typedef boost::property_tree::ptree Properties;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.main"));

class Contiging {
public:
    int run(const Properties& options) {
        if (checkOptions(options) != 0) {
            LOG4CXX_ERROR(logger, "invalid options");
            return 1;
        }
        return 0;
    }
private:
    int checkOptions(const Properties& options) {
        return 0;
    }
};

int main(int argc, char* argv[]) {
    Properties options;
    {
        // command line options
        Properties cmd;
        std::string opt_string("s:d:K:E:p:c:h");
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
