#ifndef constant_h_
#define constant_h_

#include <cstdio>

// common
extern const char* kLogConfig;
extern const char* kWorkDir;
extern const size_t kKmerSize;

// gap_filling
extern const int MAX_MERGE_GAP;
extern const int MIN_EDGE_COUNT_FOR_TRAINING;
extern const int MAX_CHOICE;
extern const int CANDI_THRESHOLD;

#endif // constant_h_
