#ifndef _kmer_hash_md5
#define _kmer_hash_md5

#include "Kmer.h"
#include "Hash_Md5.h"
#include <vector>
#include <list>
#include <deque>
#include <cmath>

using namespace std;


typedef struct Kmer_Node
{
	Kmer_Node(const Kmer &km):kmer(km){
		num = 1;
		next = NULL;
	}

	unsigned long get_code(){
		return kmer.get_code();
	}
	void show(){
		cout << kmer.get_seq() << "\t" << num << endl;
	}
	unsigned int num;
	Kmer_Node *next;
	Kmer kmer;
} Kmer_Node;
typedef Kmer_Node* Kmer_Node_Ptr;

typedef struct KNode
{
	KNode(const Kmer &km, int n=1):kmer(km),num(n){}
	KNode(const Kmer_Node &kn){
		kmer = kn.kmer;
		num = kn.num;
	}
	unsigned long get_code(){
		return kmer.get_code();
	}
	void show(){
		cout << kmer.get_seq() << "\t" << num << endl;
	}
	Kmer kmer;
	// Kmer_Node *next;
	unsigned int num;
} KNode;

class Kmer_Hash_Md5
{
public:
	Kmer_Hash_Md5(unsigned int half_first_alloc = 10240000);
	// void push_a_kmer(const Kmer&);
	// void push_a_kmer(unsigned long code);
	void push_a_kmer(unsigned long first, unsigned long second, unsigned long third);
	// void push_a_kmer(const Kmer &kmer);
	Kmer_Node_Ptr has_kmer(const Kmer &);
	// Kmer_Node_Ptr* has_kmer(const Kmer &kmer,bool &flag);
	unsigned int get_kmer_num(const Kmer&);
	unsigned long size();
	unsigned long hash_table_length();
	Kmer_Node_Ptr& operator[](unsigned long i);
	void show_all_kmers();
	void cal_lamada(double&,double&);
	void clear();
	void filter();
	void filter2();	// change buffer size
	KNode* next_node(); // traverse all node
	~Kmer_Hash_Md5();

	void test();
private:
	void rehash();
	Kmer_Hash_Md5(const Kmer_Hash_Md5 & kh);
	Kmer_Hash_Md5 & operator=(const Kmer_Hash_Md5 & kh);
private:
	static int bufferid;
	static int nodeid;
	static KNode *ret;
	unsigned int last_alloc_num;
	unsigned int last_alloc_used;

	int prime_pos;
	unsigned long hash_table_used;
	unsigned long kmer_count;
	double load_factor; 
	Kmer_Node_Ptr *hash_table;
	vector<Kmer_Node_Ptr> buffer_node;
	vector<int> buffer_size;
};
#endif
