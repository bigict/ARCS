#ifndef runner_h_
#define runner_h_

#include "config.h"

#include <iostream>
#include <map>
#include <string>
#include <tr1/tuple>
#include <vector>

#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>

typedef boost::property_tree::ptree Properties;
typedef std::vector< std::string > Arguments;

class Runner {
public:
    const std::string& options() const {
        return _options;
    }
    std::string transform(char key) const {
        std::map< char, std::string >::const_iterator i = _transform.find(key);
        if (i != _transform.end()) {
            return i->second;
        }
        return std::string(1, key);
    }
    virtual int run(const Properties& options, const Arguments& arguments) = 0;
protected:
    Runner(const std::string& options = "", const std::map< char, std::string >& table=std::map< char, std::string >()) : _options(options), _transform(table) {
    }
    std::string _options;
    std::map< char, std::string > _transform;
};

typedef Runner* RunnerPtr;

class RunnerManager {
public:
    RunnerManager() {
    }
    virtual ~RunnerManager() {
    }

    static RunnerManager* get() {
        static RunnerManager mgr;
        return &mgr;
    }

    RunnerPtr create(const std::string& name) const {
        RunnerList::const_iterator i = _runners.find(name);
        if (i != _runners.end()) {
            return std::tr1::get< 0 >(i->second);
        }
        return RunnerPtr();
    }
    bool install(const std::string& name, RunnerPtr runner, const std::string& introduction) {
        if (_runners.find(name) != _runners.end()) {
            return false;
        }
        _runners[name] = std::tr1::make_tuple(runner, introduction);
        return true;
    }
    bool uninstall(const std::string& name) {
        RunnerList::iterator i = _runners.find(name);
        if (i != _runners.end()) {
            _runners.erase(i);
            return true;
        }
        return false;
    }

    int help(int argc, char* argv[]) const {
        static std::string options("vh");
        int opt = -1;
        while ((opt = getopt(argc, argv, options.c_str())) != -1) {
            switch ((char)opt) {
            case 'v':
                std::cout << boost::format("%s version %s") % PACKAGE_NAME % PACKAGE_VERSION << std::endl;
                return 256;
            }
        }

        std::cout << boost::format("%s version %s, report bugs to %s") % PACKAGE_NAME % PACKAGE_VERSION % PACKAGE_BUGREPORT << std::endl;
        std::cout << boost::format("usage: %s <command> [<args>]") % PACKAGE_NAME << std::endl;
        std::cout << std::endl;
        std::cout << boost::format("The most commonly used %s commands are:") % PACKAGE_NAME << std::endl;

        {
            size_t max_name_length = 2;
            {
                std::vector< size_t > name_length;
                for (RunnerList::const_iterator i = _runners.begin(); i != _runners.end(); ++i) {
                    name_length.push_back(i->first.length());
                }
                if (!name_length.empty()) {
                    max_name_length += *std::max_element(name_length.begin(), name_length.end());
                }
            }
            for (RunnerList::const_iterator i = _runners.begin(); i != _runners.end(); ++i) {
                std::string cmd(i->first);
                cmd.resize(max_name_length, ' ');
                std::cout << boost::format("   %s%s") % cmd % std::tr1::get< 1 >(i->second) << std::endl;
            }
        }

        std::cout << std::endl;
        std::cout << boost::format("See '%s <command> -h' to read about a specific subcommand.") % PACKAGE_NAME << std::endl;
        return 256;
    }
private:
    typedef std::tr1::tuple< RunnerPtr, std::string > RunnerInfo;
    typedef std::map< std::string, RunnerInfo > RunnerList;
    RunnerList _runners;
};

#define RUNNER_INSTALL(name, runner, introduction) \
    RunnerManager::get()->install(name, runner, introduction)
#define RUNNER_UNINSTALL(name) \
    RunnerManager::get()->uninstall(name)

#endif // runner_h_
