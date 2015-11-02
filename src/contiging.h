#ifndef contiging_h_
#define contiging_h_

#include "runner.h"

class Contiging : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments);

private:
    Contiging();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static Contiging _runner;
};

#endif // contiging_h_
