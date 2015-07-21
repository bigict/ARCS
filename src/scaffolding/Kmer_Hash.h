#ifndef _kmer_hash
#define _kmer_hash

#include "Kmer.h"
#include <vector>
#include <list>
#include <utility>

using namespace std;


typedef struct Kmer_Node
{
	Kmer_Node():kmer_num(0){}
	Kmer kmer_buff[100000];
	unsigned int kmer_num;
} Kmer_Node;

typedef struct Kmer_Buff
{
public:
	Kmer_Buff():kmer_num(0){}
	list<Kmer_Node> content;
	unsigned long kmer_num;
} Kmer_Buff;

typedef struct Posi_Node
{
	Posi_Node():pos_num(0){}
	pair<long, long> pos_buff[100000];
	unsigned int pos_num;
} Posi_Node;

typedef struct Posi_Buff
{
public:
	Posi_Buff():pos_num(0){}
	list<Posi_Node> content;
	unsigned long pos_num;
} Posi_Buff;


class Kmer_Hash
{
public:
	void push_a_kmer(Kmer, long, long);
	void free_buff();
	void initialize_kmer_array();
	pair<long, long> get_kmer_pos(Kmer);

	unsigned long size();
	void clear();
	//	Kmer& operator[] (unsigned long pos);
	//	Kmer& at(unsigned long pos);
	
private:
	Kmer_Node kmer_node;
	Kmer_Buff kmer_buff;

	Posi_Node pos_node;
	Posi_Buff pos_buff;

	vector<Kmer> kmer_array;
	vector< pair<long, long> > pos_array;

	vector<long> kmer_index;
};
#endif