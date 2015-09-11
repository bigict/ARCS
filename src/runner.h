#ifndef runner_h_
#define runner_h_

#include <map>
#include <string>
#include <tr1/memory>

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
            return i->second;
        }
        return RunnerPtr();
    }
    bool install(const std::string& name, RunnerPtr runner) {
        if (_runners.find(name) != _runners.end()) {
            return false;
        }
        _runners[name] = runner;
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

private:
    typedef std::map< std::string, RunnerPtr > RunnerList;
    RunnerList _runners;
};

#define RUNNER_INSTALL(name, runner) \
    RunnerManager::get()->install(name, runner)
#define RUNNER_UNINSTALL(name) \
    RunnerManager::get()->uninstall(name)

#endif // runner_h_
