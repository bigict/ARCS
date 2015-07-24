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
    // command line options
    
    // config log4cxx.
    log4cxx::BasicConfigurator::configure();

    // load options.
    Properties options;
    try {
        boost::property_tree::read_ini("test.cfg", options);
    } catch (const boost::property_tree::ini_parser_error& e) {
        LOG4CXX_ERROR(logger, boost::format("load test.cfg failed(%s).") % e.what());
        return 1;
    }

    // build contigs.
    Contiging c;
    return c.run(options);
}
