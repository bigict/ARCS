#ifndef _NODE_
#define _NODE_

#include <cstring>
#include "Kmer_Short_Hash_Md5.h"
typedef struct Node
{
	Node(const Kmer_Short_Node &sn){
		kmer_short = sn.kmer_short;
		memcpy(cov, sn.num, sizeof(unsigned short)*4);
		memcpy(next, sn.nex, sizeof(unsigned int)*4);
		// cov[0] = sn.num[0];
		// cov[1] = sn.num[1];
		// cov[2] = sn.num[2];
		// cov[3] = sn.num[3];
		// next[0] = sn.nex[0];
		// next[1] = sn.nex[1];
		// next[2] = sn.nex[2];
		// next[3] = sn.nex[3];
	}
	Kmer_Short kmer_short;
	unsigned short cov[4];
	unsigned int next[4];
} Node;

#endif