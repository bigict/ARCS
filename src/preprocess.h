#ifndef preprocess_h_
#define preprocess_h_

#include "runner.h"

class Preprocess : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments);

private:
    Preprocess();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static Preprocess _runner;
};

#endif // preprocess_h_
