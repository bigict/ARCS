#include <string>
#include <vector>
#include <list>

using namespace std;

struct Dis_Node {
	size_t id;
	int dis;
	int c;
	int score;
};

struct Edge_Seq_Element {
	size_t id;
	int pos;
	int len;
};

struct Ed {
	int from;
	int to;
	int score;
};


class UniqEdgeGraph {
public:
    UniqEdgeGraph(size_t K, size_t edge_num, size_t max_overlap, size_t iteration) : _K(K), _edge_num(edge_num), _max_overlap(max_overlap), _iteration(iteration) {
    }
    virtual ~UniqEdgeGraph() {
    }

	void input_edge_pos();
	void input_edge_len();
	void input_edge_link(string);
	void input_inner_component();

	void add_a_dis(int, int, int, int,int);
	bool is_a_dis(int, int);
	void del_a_dis(int, int);
	void del_a_edge(int);
	void add_a_arc(int, int, int, int, int);
	void del_a_arc(int, int);
	void add_a_rev_arc(int, int, int, int, int);
	void del_a_rev_arc(int, int);
	int  get_dis(int, int);
	//void remove_amb_edges();
	void linearize();
	void resize(int);
	void output_graph(string);
	
	void add_a_len(int);
	void add_a_pos(int);
	Ed get_min_ed(int);
	//void add_a_ref_pos(int);
	//void choose_one();
	//void choose_forward_one(int);
	//void choose_backward_one(int);
	void initialize_component(string);
	void remove_arc_con_edge_from_overlap_pair();
	int  get_ancestor(int, int);
	int  get_descendant(int, int);
	//const vector<vector<Edge_Seq_Element> > &get_component();
	//double compute_score(int, int, int, int);
	//double compute_lam(int, int, int);
	void tran_to_line();
	void initialize_gaussian();
	double get_gaussian(int);
	int get_c(int , int);
	void output_para();

	void output_lp();

private:
	void transform_bidirection();
	void transform_single_direction();

	vector<list<Dis_Node> > arc;
	vector<list<Dis_Node> > rev_arc;
	vector< list<Dis_Node> > con;
	vector<int> pos;
	
//	vector<int> ref_pos;
	vector< vector<unsigned int> > inner_component;
	vector< vector<int> > gap;
	vector<vector<Edge_Seq_Element> > component;
	vector<vector<Edge_Seq_Element> > line_component;
	vector<int> len;
//	vector<double> gaussian;
	
	vector<pair<int, int> > overlap_pair;
	vector<int> overlap_com_id;
	//int compute_log_chimeric(int, int, int);
    size_t _K;
    size_t _edge_num;
    size_t _max_overlap;
    size_t _iteration;
};
