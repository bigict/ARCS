#ifndef debruijn_edge_h_
#define debruijn_edge_h_

#include <cassert>

#include <iostream>
#include <vector>
#include <string>

using namespace std;
 
extern int K;
extern int overlap;

// Contig
class Edge {
public:
	Edge() {
    }
	Edge(int index, long length, const std::string& seq, int copy_num) : index(index), len(length), edge_seq(seq), copy_num(copy_num) {
    }
	virtual ~Edge() {
    }

	long len;
	int copy_num;
private:
    friend std::ostream& operator<<(std::ostream &out, const Edge &obj);

	int index;	
	string edge_seq;

public:
	
	vector<int> nexts;
	vector<int> prevs;

	inline string get_edge_seq() const;

	inline string get_prefix() const;
	inline string get_suffix() const;
	
	inline string get_specific_len_prefix(  int ) const;
	inline string get_specific_len_suffix(  int ) const;

	inline char get_suf_base() const;
	inline char get_pre_base() const;
};

inline string Edge::get_edge_seq() const
{
	return edge_seq;
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
inline char Edge::get_pre_base() const {
	assert(len >= K);
	return edge_seq[K-1];
}

#endif // debruijn_edge_h_

