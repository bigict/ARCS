#ifndef contiging_h_
#define contiging_h_

#include "runner.h"

class Contiging : public Runner {
public:
    const std::string& options() const {
        return _options;
    }
    int run(const Properties& options);

private:
    Contiging() : _options("c:s:K:i:d:h") {
        RUNNER_INSTALL("contiging", this);
    }
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static Contiging _runner;
    std::string _options;
};

#endif // contiging_h_
