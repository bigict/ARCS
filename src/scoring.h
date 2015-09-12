#ifndef scoring_h_
#define scoring_h_

#include "runner.h"

class Scoring : public Runner {
public:
    int run(const Properties& options);

private:
    Scoring();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static Scoring _runner;
};

#endif // scoring_h_
