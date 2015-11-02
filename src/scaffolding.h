#ifndef scaffolding_h_
#define scaffolding_h_

#include "runner.h"

class Scaffolding : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments);

private:
    Scaffolding();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static Scaffolding _runner;
};

#endif // scaffolding_h_
