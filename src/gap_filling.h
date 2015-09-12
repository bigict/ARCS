#ifndef gap_filling_h_
#define gap_filling_h_

#include "runner.h"

class GapFilling : public Runner {
public:
    int run(const Properties& options);

private:
    GapFilling();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static GapFilling _runner;
};

#endif // gap_filling_h_
