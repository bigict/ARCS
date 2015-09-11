#ifndef runner_h_
#define runner_h_

#include "config.h"

#include <map>
#include <string>
#include <tr1/memory>
#include <tr1/tuple>

#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>

typedef boost::property_tree::ptree Properties;

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
    virtual int run(const Properties& options) = 0;
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

    int help() const {
        std::cout << boost::format("%s version %s, bugreport %s") % PACKAGE_NAME % PACKAGE_VERSION % PACKAGE_BUGREPORT << std::endl;
        std::cout << boost::format("usage: %s <command> [<args>]") % PACKAGE_NAME << std::endl;
        std::cout << std::endl;
        std::cout << boost::format("The most commonly used %s commands are:") % PACKAGE_NAME << std::endl;
        for (RunnerList::const_iterator i = _runners.begin(); i != _runners.end(); ++i) {
            std::cout << boost::format("   %s\t%s") % i->first % std::tr1::get< 1 >(i->second) << std::endl;
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
