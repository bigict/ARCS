#ifndef _hash_md5
#define _hash_md5
#define ONEG	52428800  // 100M


static const int _primes_num = 29;
static const unsigned long _prime_list[_primes_num] = {
											389,		769,
	1543,		3079,		6151,			12289,		24593,
	49157,		98317,		196613,			393241,		786433,
	1572869,	3145739,	6291469,		12582917,	25165843,
	50331653,	100663319,	201326611,		402653189,	805306457,
	1610612741,	3221225473ul,	4294967291ul, 9999999929ul, 19999999867ul,
	39999999733ul, 79999999463ul
};

#endif
