#ifndef remove_repeats_h_
#define remove_repeats_h_

#include "runner.h"

class RepeatRemover : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments);

private:
    RepeatRemover();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static RepeatRemover _runner;
};

#endif // remove_repeats_h_
