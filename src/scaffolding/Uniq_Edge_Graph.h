#ifndef uniq_edge_graph_h_
#define uniq_edge_graph_h_

#include <string>
#include <vector>
#include <list>

using namespace std;

typedef struct Dis_Node
{
	unsigned int id;
	int dis;
	int c;
	int score;
} Dis_Node;

typedef struct Edge_Seq_Element
{
	unsigned int id;
	int pos;
	int len;
} Edge_Seq_Element;

class Uniq_Edge_Graph
{
public:
	void add_a_dis(int, int, int, int);
	bool is_a_dis(int, int);
	void del_a_dis(int, int);
	int get_dis(int, int);
	void linearize();
	void resize(int);
	void output_graph(string);
	void output_detailed_graph(string);	
	double compute_score(int, int, int, int);
	void initialize_gaussian();
	int get_c(int , int);
	void output_lp();
	void output_edge_len();
	void set_edge_score();
	void set_prob_cutoff();
	void remove_amb_edges();

	/*  Zheng quangang
	 * date : 2015-03-30
	 */
	void build_long_contig_graph();
	unsigned int next_contig(unsigned int);
	void remove_links(unsigned int);
	void extend_lcg();
	void input_real_dis();
	/*
	 * End
	 */

private:
	int log_pro_threshold;
};

#endif // uniq_edge_graph_h_
