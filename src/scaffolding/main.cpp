#include "linearization.h"

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

typedef boost::property_tree::ptree Properties;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffolding.main"));

using namespace std;

extern int K;
extern int PAIR_KMER_CUTOFF;
extern int PAIR_READS_CUTOFF;
extern int EDGE_CUTOFF;
extern int PAIR_KMER_NUM ;
extern int GENOME_LEN;
extern int INSERT_SIZE ;
extern double DELTA ;
extern int CPU_NUM;
extern double LINK_QUALITY_PERCENT ;
extern int LOG_PRO_THRESHOLD ;
extern int ITER;

int K = 31;
int PAIR_KMER_CUTOFF = 10;
int PAIR_READS_CUTOFF = 5;
int EDGE_CUTOFF = 0;
int PAIR_KMER_NUM = 0;
int GENOME_LEN = 4000000;
int INSERT_SIZE = 180;
double DELTA = 5.0;
int CPU_NUM = 8;
double LINK_QUALITY_PERCENT = 0.01;
int LOG_PRO_THRESHOLD = -200;
int ITER = 0;

static const char *optString = "d:K:c:e:L:D:1:2:P:p:i:r:R:";

class Scaffolding {
public:
    int run(const Properties& options) {
        if (checkOptions(options) != 0) {
            return 1;
        }

        Linearization li(options.get< size_t >("K"));
        li.initialize_edge(options.get< std::string >("c"));
        li.initialize_scaf();
        li.linearize(options.get< std::string >("1"), options.get< std::string >("2"));

        return 0;
    }

private:
    int checkOptions(const Properties& options);
    int printHelps() const;
};

int Scaffolding::checkOptions(const Properties& options) {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }

    std::vector< std::string > necessary = boost::assign::list_of("c")("1")("2")("K");
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

int Scaffolding::printHelps() const {
    return 1;
}

int main(int argc, char **argv) {
    Properties options;
    {
        // command line options
        std::string opt_string("d:K:c:e:L:D:1:2:P:p:i:r:R:h");
        int opt = -1;
        while ((opt = getopt(argc, argv, opt_string.c_str())) != -1) {
            std::string key(1, (char)opt);
            if (optarg == NULL) {
                options.put(key, NULL);
            } else {
                options.put(key, optarg);
            }
        }
    }
	
    Scaffolding s;
    return s.run(options);
}

