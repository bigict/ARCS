#ifndef _kmer_hash
#define _kmer_hash

#include "Kmer.h"
#include <vector>
#include <list>

using namespace std;

typedef struct Buff_Node
{
	Buff_Node():kmer_num(0){}
	Kmer kmer_buff[100000];
	unsigned int kmer_num;
} Buff_Node;

class Kmer_Buff
{
public:
	Kmer_Buff():kmer_num(0){}
	list<Buff_Node> content;
	unsigned long kmer_num;
};

class Kmer_Hash
{
public:
	void push_a_kmer(Kmer);
	void free_buff();
	void initialize_kmer_array();
	unsigned int get_kmer_num(Kmer);
	unsigned int get_kmer_num(unsigned long);
	unsigned long size();
	void clear();
	Kmer& operator[] (unsigned long pos);
	Kmer& at(unsigned long pos);
	void filter();
private:
	Buff_Node buff_node;
	Kmer_Buff kmer_buff;
	vector<Kmer> kmer_array;
	vector<unsigned int> kmer_num;
};
#endif
