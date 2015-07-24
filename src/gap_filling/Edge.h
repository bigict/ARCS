#ifndef _EDGE_H
#define _EDGE_H

#include <iostream>
#include <vector>
#include <string>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


#include "Constant.h"

using namespace std;
 
extern int K;
extern int overlap;

//Contig
class Edge
{

friend ostream &operator<<(ostream &out, const Edge &obj);

private:

	static const int C_IN_EDGE = 4;

	int index;	
	long len;
	int copy_number;
	string edge_seq;

public:
	
	vector<int> nexts;
	vector<int> prevs;

	Edge();
	Edge(int, long, string, int);

	~Edge();

	inline long get_len() const;
	inline string get_edge_seq() const;
	inline int get_copy_number() const;
	inline void sub_copy_number() ;

	inline string get_prefix() const;
	inline string get_suffix() const;
	
	inline string get_specific_len_prefix(  int ) const;
	inline string get_specific_len_suffix(  int ) const;

	inline char get_suf_base() const;
	inline char get_pre_base() const;
};

inline long Edge::get_len() const
{
	return len;
}

inline string Edge::get_edge_seq() const
{
	return edge_seq;
}

inline int Edge::get_copy_number() const
{
	return copy_number;
}

inline void Edge::sub_copy_number() 
{
	--copy_number;
}

inline string Edge::get_prefix() const
{
	assert(len >= K);
	string pre_kmer = edge_seq.substr(0,K-1);
	return pre_kmer;
}

inline string Edge::get_suffix() const
{
	assert(len >= K);
	string suf_kmer = edge_seq.substr(len -(K-1), K-1);
	return suf_kmer;
}

inline string Edge::get_specific_len_prefix(  int overlap) const
{
	assert(overlap >= 0);
	string pre_kmer = edge_seq.substr(0,overlap);
	return pre_kmer;
}

inline string Edge::get_specific_len_suffix(  int overlap) const
{
	assert(overlap >= 0);
	string suf_kmer = edge_seq.substr( len - overlap, overlap);
	return suf_kmer;
}

inline char Edge::get_suf_base() const
{
	assert(len >= K);
	return edge_seq[len-K];
}
inline char Edge::get_pre_base() const
{
	assert(len >= K);
	return edge_seq[K-1];
}

#endif

