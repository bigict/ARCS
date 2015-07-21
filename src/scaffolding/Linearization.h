#include <string>

#include "Pos_Recorder.h"
#include "Edge.h"
#include "Kmer_Hash.h"
#include "Lib.h"

using namespace std;

extern int K;

class Linearization
{
public:
	void initialize_kmer_hash();
	void initialize_kmer_hash_for_trainning();
	void free_kmer_hash();
	void initialize_edge(string);
	void initialize_scaf();
	void map_read(string, string);
	void linearize(string, string);
	void print_edge();
	void output_para();
	void estimate_insert_size();
private:
	Uniq_Edge_Graph u_e_g;
	
	vector<Edge> edge;	

	vector< vector<unsigned int> > component;
	vector< vector<unsigned int> > gap;
};

