#ifndef preprocess_h_
#define preprocess_h_

#include "runner.h"

class Preprocess : public Runner {
public:
    const std::string& options() const {
        return _options;
    }
    int run(const Properties& options);

private:
    Preprocess() : _options("c:s:K:n:i:o:h") {
        RUNNER_INSTALL("preprocess", this);
    }
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static Preprocess _runner;
    std::string _options;
};

#endif // preprocess_h_
