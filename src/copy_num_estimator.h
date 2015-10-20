#ifndef copy_num_estimator_h_
#define copy_num_estimator_h_

#include "runner.h"

class CopyNumEstimator : public Runner {
public:
    int run(const Properties& options, const Arguments& arguments);

private:
    CopyNumEstimator();
    int checkOptions(const Properties& options) const;
    int printHelps() const;

    static CopyNumEstimator _runner;
};

#endif // copy_num_estimator_h_
