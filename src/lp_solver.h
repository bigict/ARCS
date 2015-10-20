#ifndef lp_solver_h_
#define lp_solver_h_

#include "runner.h"

class LPSolver : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments);

private:
    LPSolver();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static LPSolver _runner;
};

#endif // lp_solver_h_
