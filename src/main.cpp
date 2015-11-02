#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

#include "constant.h"
#include "runner.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.main"));

int main(int argc, char* argv[]) {
    if (argc < 2) {
        return RunnerManager::get()->help(argc, argv);
    }

    RunnerPtr runner = RunnerManager::get()->create(argv[1]);
    if (!runner) {
        return RunnerManager::get()->help(argc, argv);
    }

    Properties options;
    {
        // command line options
        Properties cmd;
        const std::string& opt_string = runner->options();
        int opt = -1;
        while ((opt = getopt(argc - 1, argv + 1, opt_string.c_str())) != -1) {
            std::string key = runner->transform((char)opt);
            if (optarg == NULL) {
                cmd.put(key, NULL);
            } else {
                std::string val = optarg;
                if (cmd.find(key) != cmd.not_found()) {
                    val = boost::str(boost::format("%s:%s") % cmd.get< std::string >(key) % val);
                }
                cmd.put(key, val);
            }
        }

        // config log4cxx.
        const std::string log_config = cmd.get< std::string >("c", kLogConfig);
        if (boost::filesystem::exists(log_config)) {
            log4cxx::PropertyConfigurator::configure(log_config);
        } else {
            log4cxx::BasicConfigurator::configure();
        }
        
        // load ini options
        if (cmd.find("s") != cmd.not_found()) {
            const std::string file_config = cmd.get< std::string >("s");
            try {
                boost::property_tree::read_ini(file_config, options);
            } catch (const boost::property_tree::ini_parser_error& e) {
                LOG4CXX_ERROR(logger, boost::format("load %s failed(%s).") % file_config % e.what());
                return 1;
            }
        }

        // merge options
        for (Properties::const_iterator it = cmd.begin(); it != cmd.end(); it++){
            options.put(it->first,it->second.data());
        }
    }

    Arguments arguments;
    for (int i = optind + 1; i < argc; ++i) {
        arguments.push_back(argv[i]);
    }
    //std::copy(argv, optind, argc, std::back_inserter(arguments));

    return runner->run(options, arguments);
}
