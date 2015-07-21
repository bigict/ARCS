#include <string>
#include <pthread.h>
#include "Kmer_Hash_Md5.h"
#include "Lib.h"
#include "Kmer_Short_Hash_Md5.h"
#include "Kmer.h"
#include "Kmer_Short.h"
#include "node.h"

using namespace std;

extern int K;
const int BUFFER_SIZE = 60000;

class DeBruijnGraph;

typedef struct ARG{
	pthread_t threadi;
	int threshold;
	int minthres;
	float percent;
	ifstream *left;
	ifstream *right;
	unsigned long (*buffer)[BUFFER_SIZE][3];
	pthread_mutex_t *buffer1;
	pthread_mutex_t *buffer2;
	pthread_cond_t *cond1;
	pthread_cond_t *cond2;
	DeBruijnGraph *db;
	Kmer_Hash_Md5 *kmer_seq;
}ARG;

// typedef struct Node
// {
// 	Kmer_Short kmer_short;
// 	int cov[4];
// 	int next[4];
// } Node;

typedef struct Long_Node
{
	vector<char> edge_seq;
	vector<int> cov;
	int next;
//	Long_Node *next_node;
} Long_Node;

class DeBruijnGraph
{
public:
//	DeBruijnGraph();
	void initialize_kmer(Lib &lib);
	void output_kmer();

	void initialize_graph();
	void output_graph(string tmp);
	
	void build_condensed_de_bruijn_graph();
	void split_identical_edge();
	void output_min_cost_flow(string file_name);
	string reverse_complement(const string &fw);
	string reverse_complement_replace_N(string &fw);

	void add_con_edge(int index, vector<char> &c_edge_seq, vector<int> &c_cov, int c_next);
	void output_con_graph();

	void trim_bubble_tip();
	void delete_con_edge(list<int> &node_index, list<int> &node_next_index);

	void set_lamada();
	void set_n();
	void output_parameter(string);

	void output_condensed_de_Bruijn_graph(string);
	void output_de_Bruijn_graph(string file_name);

private:
	int kmer_size;

	Kmer_Hash_Md5 kmer_seq;
	Kmer_Short_Hash_Md5 kmer_short;
	
	vector<Node> graph;
	vector<vector<Long_Node> > con_graph;

	double lamada;
	double n;
	double N;
private:
	double get_cost(vector<int> &cov, double di);
	double get_d0(vector<int> &cov);
	double get_k1(vector<int> &cov);
	double get_k2(vector<int> &cov);

	void divide_edge(int, int);
};
