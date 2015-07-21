#ifndef _kmer_short_hash_md5
#define _kmer_short_hash_md5

#include "Kmer_Short.h"
#include "Hash_Md5.h"
#include <vector>
#include <list>
#include <cmath>

using namespace std;

typedef struct Kmer_Short_Node
{
	Kmer_Short_Node(const Kmer_Short &km, int index=0, unsigned int numb=0):kmer_short(km){
		num[0] = 0;num[1] = 0;num[2] = 0;num[3] = 0;
		num[index] = numb;
		pos = 1;
		next = NULL;
	}

	void set(int index, unsigned int numb){
		num[index] = numb;
	}

	unsigned long get_code(){
		return kmer_short.get_code();
	}
	void show(){
		cout << kmer_short.get_seq() << "\t" << num << endl;
	}
	unsigned short num[4];
	long pos;
	unsigned int nex[4];
	Kmer_Short kmer_short;
	Kmer_Short_Node *next;
} Kmer_Short_Node;
typedef Kmer_Short_Node* Kmer_Short_Node_Ptr;

#include "node.h"

class Kmer_Short_Hash_Md5
{
public:
	Kmer_Short_Hash_Md5(unsigned int half_first_alloc=512);
	void push_a_kmer_short(const Kmer_Short&);
	void push_a_kmer_short(const Kmer&, unsigned int numb);
	Kmer_Short_Node_Ptr has_kmer(const Kmer_Short &);
	unsigned short get_kmer_num(const Kmer_Short&,int);
	unsigned short* get_kmer_nums(const Kmer_Short&);
	unsigned long size();
	unsigned long hash_table_length();
	Kmer_Short_Node_Ptr& operator[](unsigned long i);
	void show_all_kmers();
	void clear();
	void filter();
	void esti_pos();
	void transform2graph(vector<Node>&);
	~Kmer_Short_Hash_Md5();
private:
	void rehash();
	Kmer_Short_Hash_Md5(const Kmer_Short_Hash_Md5 & kh);
	Kmer_Short_Hash_Md5 & operator=(const Kmer_Short_Hash_Md5 & kh);
private:
	int buffers;
	vector<Kmer_Short_Node_Ptr> buffer_node;
	vector<int> buffer_size;
	unsigned int last_alloc_num;
	unsigned int last_alloc_used;

	Kmer_Short_Node_Ptr *hash_table;
	int prime_pos;
	unsigned long hash_table_used;
	unsigned long kmer_count;
	double load_factor; 
};
#endif
