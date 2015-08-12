#include <map>
#include <vector>
#include <string>
#include <list>
#include <pthread.h>

#include "uniq_edge_graph.h"

using namespace std;

typedef struct Pos_Node
{
	int id;
	long dis_sum;
	int c;
	int reads;
} Pos_Node;

class Pos_Recorder
{
public:
	void add_a_dis(int, int, int, int,int);
	void compute_dis(Uniq_Edge_Graph &);
	void resize(int);
	void output_pos_graph();
	void initialize_mutex_array(int);
	void destroy_mutex_array();
private:
	vector< list<Pos_Node> > con;
	vector<pthread_mutex_t> mutex_array;
};
