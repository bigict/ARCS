#ifndef lp_divider_h_
#define lp_divider_h_

#include "runner.h"

class LPDivider : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments);

private:
    LPDivider();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static LPDivider _runner;
};

#endif // lp_divider_h_
