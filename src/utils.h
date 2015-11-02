#ifndef utils_h_
#define utils_h_

#ifndef SIZEOF_BITS
#define SIZEOF_BITS(x) (8 * sizeof(x))
#endif

#ifndef SIZEOF_ARRAY
#define SIZEOF_ARRAY(x)  (sizeof(x) / sizeof(x[0]))
#endif

#endif // utils_h_
