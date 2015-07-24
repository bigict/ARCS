#ifndef _kmer_short_hash
#define _kmer_short_hash

#include"Kmer_Short.h"

#include <vector>
#include <list>

using namespace std;

typedef struct Short_Buff_Node
{
	Short_Buff_Node():kmer_short_num(0){}
	Kmer_Short kmer_short_buff[100000];
	unsigned int kmer_short_num;
} Short_Buff_Node;

class Kmer_Short_Buff
{
public:
	Kmer_Short_Buff():kmer_short_num( 0){}
	list<Short_Buff_Node> content;
	unsigned long kmer_short_num;
};

class Kmer_Short_Hash
{
public:
	void push_a_kmer_short(Kmer_Short);
	void free_buff();
	void initialize_kmer_short_array();
	long get_kmer_short_pos(Kmer_Short);
	unsigned long size();
	void clear();
	Kmer_Short& operator[] (unsigned long pos);
	Kmer_Short& at(unsigned long pos);
private:
	Short_Buff_Node short_buff_node;
	Kmer_Short_Buff kmer_short_buff;
	vector<Kmer_Short> kmer_short_array;
};

#endif
