/***************************************************************************
 * Description:
 *          *
 *
 string DeBruijnGraph::reverse_complement(const string &fw)
 *
 void DeBruijnGraph::add_con_edge(int index, vector<char> &c_edge_seq, vector<int> &c_cov, int c_next)
 *
 void DeBruijnGraph::deal_identical_edge()
 *
 void DeBruijnGraph::convert()
 *
 void DeBruijnGraph::divide_edge(int i, int j)
 *
 * Author  : Renyu Zhang
 * Language: C
 * Date    : 2014-04-01
 ***************************************************************************/

#include <stdlib.h>
#include <list>
#include <fstream>
#include <assert.h>
#include <cmath>
#include <pthread.h>
#include <algorithm>
#include "DeBruijnGraph.h"

extern int NUM_THREAD;
extern int K;
extern bool FILTER;
extern string rev_tem("TGAC");
extern int READ_LENGTH_CUTOFF;

string DeBruijnGraph::reverse_complement(const string &fw)
{
	string rc;
	rc.resize(fw.size());
	int len = fw.size();
	for (int i = len - 1; i >= 0; i--)
	{
		rc[len - 1 - i] = rev_tem[(fw.at(i) & 0x06) >> 1];
	}
	return rc;
}
string DeBruijnGraph::reverse_complement_replace_N(string &fw)
{
	string rc;
	rc.resize(fw.size());
	int len = fw.size();
	for (int i = len - 1; i >= 0; i--)
	{
		if(fw.at(i) == 'N')
			fw[i] = rand() & 0x06;
		rc[len - 1 - i] = rev_tem[(fw.at(i) & 0x06) >> 1];
	}
	return rc;
}

void DeBruijnGraph::add_con_edge(int index, vector<char> &c_edge_seq, vector<int> &c_cov, int c_next)
{
	int max;
	if(index > c_next)
		max = index;
	else
		max = c_next;

	if (index >= con_graph.size())
		con_graph.resize(max+1);
	Long_Node new_edge;
	new_edge.edge_seq = c_edge_seq;
	new_edge.cov = c_cov;
	new_edge.next = c_next;

	con_graph[index].push_back(new_edge);

}

/***************************************************************************
 * Description:
 *          
 * Author  : wangbing, zheng quangang
 * Language: C ++
 * Date    : 2015-02-09, 2015-04-09
 ***************************************************************************/

pthread_mutex_t mutex;
pthread_mutex_t print_mutex;

vector<int> initial_node;
vector<int> new_pos;
vector<int> in_degree;
vector<int> out_degree;
unsigned int sum_edge;
struct send_to_thread{
	int i;
	DeBruijnGraph *de;
};

void * call_back(void * argc){

	DeBruijnGraph* de = ((send_to_thread*)argc)->de;
	int num_thread = ((send_to_thread*)argc)->i;
	int start=0, end=0;
	start = num_thread*(initial_node.size() / NUM_THREAD);
	end = start + (initial_node.size() / NUM_THREAD);
	if (num_thread == NUM_THREAD-1){
		end = initial_node.size();
	}
	int size = end - start;

	// pthread_mutex_lock( &print_mutex );
	// cout << "start " << num_thread << " thread for condence graph" << endl;
	// pthread_mutex_unlock( &print_mutex );

	if(size == 0)
		pthread_exit(NULL);

	vector<char> tmp_str;
	vector<int> tmp_cov;

	vector<int>::iterator it;
	int i = start;
	string tmp_kmer_seq;

	const char tem[4] = {'A', 'C', 'T', 'G'};		//"ACTG";
	bool flag = false;
	int count = 0;
	int index;
	it = initial_node.begin() + (unsigned int)start;

	for(; it != initial_node.begin()+(unsigned int)end; it++)
	{
		for(int j = 0; j < 4; j++)
		{
			index = *it;
			if (de->graph[index].cov[j] > 0)
			{
				tmp_kmer_seq = de->graph[index].kmer_short.get_seq();

//				assert(tmp_kmer_seq.size() == K - 1);

//				for(int kmer_seq_i = 0; kmer_seq_i < tmp_kmer_seq.size(); kmer_seq_i ++)
//				{
//					tmp_str.push_back( tmp_kmer_seq[kmer_seq_i]);
//				}
				const char *p = tmp_kmer_seq.c_str();
				tmp_str.insert(tmp_str.end(),p,p+tmp_kmer_seq.length());

				tmp_str.push_back(tem[j]);
				tmp_cov.push_back(de->graph[index].cov[j]);

				index = de->graph[index].next[j];

				while(in_degree[index] == 1 && out_degree[index] == 1)
				{
					//flag = false;
					for(int tmp_j = 0; tmp_j < 4; tmp_j++)
					{
						if (de->graph[index].cov[tmp_j] > 0)
						{	
							tmp_str.push_back(tem[tmp_j]);
							tmp_cov.push_back(de->graph[index].cov[tmp_j]);

							index = de->graph[index].next[tmp_j];
							//flag = true;
							break;
						}
					}
				}

				pthread_mutex_lock( &mutex );
				de->add_con_edge(i, tmp_str, tmp_cov, new_pos[index]);
				sum_edge ++;
				pthread_mutex_unlock( &mutex );

				tmp_str.clear();
				tmp_cov.clear();
			}	
		}
		i ++;
	}
	pthread_exit(NULL);
}

void DeBruijnGraph::build_condensed_de_bruijn_graph()
{
	cout << "Building condensed de bruijn graph" << endl;

	con_graph.clear();
	in_degree.clear();
	in_degree.resize(graph.size());
	out_degree.clear();
	out_degree.resize(graph.size());

	fill(in_degree.begin(), in_degree.end(), 0);
	fill(out_degree.begin(), out_degree.end(), 0);

	for(int i = 0; i < graph.size(); i++)
	{
		for(int j = 0; j < 4; j++)
		{
			if (graph[i].cov[j] > 0)
			{
				in_degree[graph[i].next[j]] ++;
				out_degree[i] ++;
			}
		}
	}
	initial_node.clear();

	new_pos.resize(graph.size());
	fill(new_pos.begin(), new_pos.end(), -1);
	int index = 0;
	for(int i = 0; i < graph.size(); i ++)
	{
		if((in_degree[i] != 1 || out_degree[i] != 1) && out_degree[i] != 0 )
		{
			initial_node.push_back(i);
			new_pos[i] = index;
			index ++;
		}
	}

	con_graph.resize(initial_node.size());
	sum_edge = 0;
	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&print_mutex, NULL);
	pthread_t thr[NUM_THREAD];
	send_to_thread p[NUM_THREAD];
	for(int i=0; i<NUM_THREAD; i++){
		p[i].i = i;
		p[i].de = this;
		int ret = 0;
		if( ret = pthread_create(&(thr[i]), NULL, call_back, (void*)&p[i]) )
		{
			cout << "[ERROR] create thread " << i << " error. error code " << ret << endl;
			exit(1);
		}
	}
	cout << "[INFO] " << NUM_THREAD << " threads for condence graph" << endl;
	for(int i=0; i<NUM_THREAD; i++){
		pthread_join(thr[i], NULL);
	}
	pthread_mutex_destroy(&mutex);
	pthread_mutex_destroy(&print_mutex);
	cout << "\tcondensed node number = " << con_graph.size() << endl;
	cout << "\tcondensed edge number = " << sum_edge << endl;
}

/* end wangbing*/
/*
void DeBruijnGraph::build_condensed_de_bruijn_graph()
{
	cout << "Building condensed de bruijn graph" << endl;

	con_graph.clear();

	vector<int> in_degree;
	in_degree.resize(graph.size());
	vector<int> out_degree;
	out_degree.resize(graph.size());

	for(int i = 0; i < graph.size(); i++)
	{
		in_degree[i] = 0;
		out_degree[i] = 0;
	}

//	cout << "\tcalculate in_degree out_degree" << endl;
	for(int i = 0; i < graph.size(); i++)
	{
		for(int j = 0; j < 4; j++)
		{
			if (graph[i].cov[j] > 0)
			{
				in_degree[graph[i].next[j]] ++;
				out_degree[i] ++;
			}
		}
	}

	list<int> initial_node;
	vector<int> new_pos;
	new_pos.resize(graph.size());

	for(int i = 0; i < graph.size(); i++)
	{
		new_pos[i] = -1;
	}

//	cout << "\tproduce initial node" << endl;
	int index = 0;
	for(int i = 0; i < graph.size(); i ++)
	{
		if((in_degree[i] != 1 || out_degree[i] != 1) && !(in_degree[i] == 0 && out_degree[i] == 0) )
		{
			initial_node.push_back(i);
			new_pos[i] = index;
			index ++;
		}
	}

	//cout << "initial node produce end!!! " << index << endl;

	con_graph.resize(initial_node.size());

	vector<char> tmp_str;
	vector<int> tmp_cov;

	const char tem[4] = {'A', 'C', 'T', 'G'}; 
	//"ACTG";

	list<int>::iterator it;
	int i = 0;
	string tmp_kmer_seq;

	bool flag = false;
	unsigned int sum = 0;
	for(it = initial_node.begin(); it != initial_node.end(); it++)
	{
		for(int j = 0; j < 4; j++)
		{

			index = *it;
			if (graph[index].cov[j] > 0)
			{
				tmp_kmer_seq = graph[index].kmer_short.get_seq();

//				assert(tmp_kmer_seq.size() == K - 1);

//				for(int kmer_seq_i = 0; kmer_seq_i < tmp_kmer_seq.size(); kmer_seq_i ++)
//				{
//					tmp_str.push_back( tmp_kmer_seq[kmer_seq_i]);
//				}
				const char *p = tmp_kmer_seq.c_str();
				tmp_str.insert(tmp_str.end(),p,p+tmp_kmer_seq.length());

				tmp_str.push_back(tem[j]);
				tmp_cov.push_back(graph[index].cov[j]);

				index = graph[index].next[j];

				while(in_degree[index] == 1 && out_degree[index] == 1)
				{
					flag = false;
					for(int tmp_j = 0; tmp_j < 4; tmp_j++)
					{
						if (graph[index].cov[tmp_j] > 0)
						{	
							tmp_str.push_back(tem[tmp_j]);
							tmp_cov.push_back(graph[index].cov[tmp_j]);

							index = graph[index].next[tmp_j];
							flag = true;
							break;
						}
					}
				}


				add_con_edge(i, tmp_str, tmp_cov, new_pos[index]);
				sum ++;

				tmp_str.clear();
				tmp_cov.clear();
			}	
		}
		i ++;
	}

	cout << "\tcondensed node number = " << con_graph.size() << endl;
	cout << "\tcondensed edge number = " << sum << endl;
	con_num = sum;
}
*/
void DeBruijnGraph::split_identical_edge()
{
	vector<vector<int> > vec_pos;
	int cur;
	int k;
	int num = 0;
	for (int i = 0; i < con_graph.size(); i++)
	{
		for (int j = 0; j < con_graph[i].size(); j++)
		{
			cur = con_graph[i][j].next;
			for ( k = 0; k < vec_pos.size(); k++)
			{
				if( con_graph[i][vec_pos[k][0]].next == cur )
				{
					break;
				}
			}
			if ( k < vec_pos.size())
			{
				vec_pos[k].push_back(j);
			}else
			{
				vector<int> tmp;
				tmp.push_back(j);
				vec_pos.push_back(tmp);
			}
		}

		for (int j = 0; j < vec_pos.size(); j++)
		{
			if (vec_pos[j].size() > 1)
			{
				for (int tmp_k = 1; tmp_k < vec_pos[j].size(); tmp_k ++)
				{
					divide_edge(i,vec_pos[j][tmp_k]);
					num ++;
				}
			}
		}
		vec_pos.clear();
	}
	cout << "[TEST]split identical edges = " << num << endl;
}

void DeBruijnGraph::divide_edge(int i, int j)
{
	vector<int> fir_cov, sec_cov;
	vector<char> fir_str, sec_str;
	int media;
	int size;
	int next;

	media = con_graph[i][j].cov.size()/2;
	if (media == 0)
		return ;

	for ( int tmp_i = 0; tmp_i < media; tmp_i++)
	{
		fir_cov.push_back(con_graph[i][j].cov[tmp_i]);
	}
	for ( int tmp_i = media; tmp_i < con_graph[i][j].cov.size(); tmp_i++)
	{
		sec_cov.push_back(con_graph[i][j].cov[tmp_i]);
	}
	for ( int tmp_i = 0; tmp_i < media + K - 1; tmp_i++)
	{
		fir_str.push_back(con_graph[i][j].edge_seq[tmp_i]);
	}
	for ( int tmp_i = media; tmp_i < con_graph[i][j].cov.size() + K - 1; tmp_i++)
	{
		sec_str.push_back(con_graph[i][j].edge_seq[tmp_i]);
	}
	next = con_graph[i][j].next;
	size = con_graph.size();

	con_graph[i][j].edge_seq.swap(fir_str);
	con_graph[i][j].cov.swap(fir_cov);
	con_graph[i][j].next = size;
	add_con_edge(size, sec_str, sec_cov, next);
}

double avg(list<int> &cov)
{
	if( cov.size() == 0)
		return 0;
	double sum = 0;
	for(list<int>::iterator it = cov.begin(); it != cov.end(); it++)
	{
		sum += *it;
	}
	return (sum/cov.size());
}


void DeBruijnGraph::delete_con_edge(list<int> &node_index, list<int> &node_next_index)
{

	if(node_index.size() != node_next_index.size())
	{
		cout << " Delete failed ... size not ok" << endl;
		return;
	}
	list<int>::iterator it_i = node_index.begin();
	list<int>::iterator it_n = node_next_index.begin();
	while(it_i != node_index.end())
	{
		graph[*it_i].cov[*it_n] = 0;
		graph[*it_i].next[*it_n] = 0;
		it_i ++;
		it_n ++;
	}

}
/*
void DeBruijnGraph::set_lamada()
{
	cout << "Calculating k-mer coverage distribution" << endl;
	int count = 0;
	double sum = 0;
	int cov;

	vector<int> cov_vec;

	for(unsigned long it = kmer_seq.size()/2 - 50000; it < kmer_seq.size(); it++ )
	{
		cov = kmer_seq.get_kmer_num(it);
		if (cov > 3.0)
		{
			sum += cov;
			count ++;
			cov_vec.push_back(cov);
			if(count > 100000)
			{
				break;
			}
		}
	}
	
	lamada = sum / count;

	if (lamada > 100)
	{
		count = 0;
		sum = 0;
		
		cov_vec.clear();

		for(unsigned long it = kmer_seq.size() /2 - 50000; it < kmer_seq.size(); it++ )
		{
			cov = kmer_seq.get_kmer_num(it);
			if (cov > lamada/20.0 && cov < 2*lamada)
			{
				sum += cov;
				count ++;
				cov_vec.push_back(cov);
				if(count > 100000)
				{
					break;
				}
			}

		}
		lamada = sum / count;
	}
	sum = 0;
	for (int i = 0; i < cov_vec.size(); i++)
	{
		sum += pow(cov_vec[i] - lamada, 2);
	}
	sum = sum/(cov_vec.size() - 1);
	double var = sqrt(sum);

	cout << "\tmean = " << lamada << endl;
	cout << "\tvariance = " << var << endl;
	
	if (var > 1.7 * lamada)
	{
		cout << "Warning: variance is too large" << endl;
	}

}
*/
/*********************
 * Author: Zheng
 * Date  : 2015-3-3
 ********************/
void DeBruijnGraph::set_lamada()
{
	cout << "Calculating k-mer coverage distribution" << endl;
	double var = 0.0;
	kmer_seq.cal_lamada(lamada,var);
	cout << "\tmean = " << lamada << endl;
	cout << "\tvariance = " << var << endl;

	if (var > 1.7 * lamada)
	{
		cout << "Warning: variance is too large" << endl;
	}
}

void DeBruijnGraph::set_n()
{
	double sum = 0;

	for(int i = 0; i < graph.size(); i++)
	{
		for(int j = 0; j < 4; j++)
		{
			if (graph[i].cov[j] > 0)
			{
				sum += graph[i].cov[j];
			}
		}
	}

	n = sum;
	N = n / lamada;
}

void DeBruijnGraph::trim_bubble_tip()
{
	cout << "Removing tips, bubbles, and some chimeric paired reads" << endl;
	vector<int> in_degree;
	in_degree.resize(graph.size());
	vector<int> out_degree;
	out_degree.resize(graph.size());
	unsigned long graph_size = graph.size();
	for(int i = 0; i < graph_size; i++)
	{
		in_degree[i] = 0;
		out_degree[i] = 0;				
	}
//	cout << "graph.size() = " << graph.size() << endl;
//	cout << "\tcalculate in_degree out_degree" << endl;
	for(int i = 0; i < graph_size; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			if (graph[i].cov[j] > 0)
			{
				in_degree[graph[i].next[j]] ++;
				out_degree[i] ++;
			}
		}
	}

	list<int> initial_node;
	vector<int> new_pos;
	new_pos.resize(graph_size);
	for(int i = 0; i < graph_size; i++)
	{
		new_pos[i] = -1;
	}

	int index = 0;
	for(int i = 0; i < graph_size; i ++)
	{
		if (in_degree[i] != 1 || out_degree[i] != 1)
		{
			initial_node.push_back(i);
			new_pos[i] = index;
			index ++;
		}
	}

	list<int> node_index;
	list<int> node_next_index;
	list<int> node_cov;

	unsigned int sum = 0;
	unsigned int sum1 = 0;
	list<int>::iterator it;
	for(it = initial_node.begin(); it != initial_node.end(); it++)
	{
		for(int j = 0; j < 4; j++)
		{
			index = *it;
			if (graph[index].cov[j] > 0)
			{
				node_index.push_back(index);
				node_next_index.push_back(j);
				node_cov.push_back(graph[index].cov[j]);
				index = graph[index].next[j];

				while(in_degree[index] == 1 && out_degree[index] == 1)
				{
					for(int tmp_j = 0; tmp_j < 4; tmp_j++)
					{
						if (graph[index].cov[tmp_j] > 0)
						{	
							node_index.push_back(index);
							node_next_index.push_back(tmp_j);
							node_cov.push_back(graph[index].cov[tmp_j]);
							index = graph[index].next[tmp_j];
							break;
						}
					}
				}

				//if(in_degree[node_index.front()] == 0 || out_degree[graph[node_index.back()].next[node_next_index.back()]] == 0 )
				if((in_degree[node_index.front()] == 0 && out_degree[node_index.front] == 1) || (in_degree[graph[node_index.back().next[node_next_index.back()]]] == 1 && out_degree[graph[node_index.back()].next[node_next_index.back()]] == 0 ))
				{
					if(node_index.size()+K-1 < 2*K)
					{
						delete_con_edge(node_index, node_next_index);
						sum ++;
					}else if(avg(node_cov) < lamada/5.0)
					{
						delete_con_edge(node_index, node_next_index);
						sum ++;
					}
				}
				else if(node_index.size()+K-1 < 2*K && avg(node_cov) < lamada/5.0)
				{
					delete_con_edge(node_index, node_next_index);
					list<int>::iterator it_i = --(node_index.end());
					list<int>::iterator it_n = --(node_next_index.end());
					graph[*it_i].cov[*it_n] = 0;
					graph[*it_i].next[*it_n] = 0;
					sum1 ++;
				}
				else if(avg(node_cov) <= 3.0 || avg(node_cov) < lamada/10.0)
				{
					list<int>::iterator it_i = --(node_index.end());
					list<int>::iterator it_n = --(node_next_index.end());
					graph[*it_i].cov[*it_n] = 0;
					graph[*it_i].next[*it_n] = 0;
					sum1 ++;
				}

				node_index.clear();
				node_next_index.clear();
				node_cov.clear();
			}//if	
		}//for j 
		//	i++;
	}//for node

	cout << "\tRemoved condensed edge number(Tips) = " << sum << endl;

	cout << "\tRemoved condensed edge number(Inner edges) = " << sum1 << endl;
}

void DeBruijnGraph::output_con_graph()
{
	cout << "output con graph" << endl;

	for(int i = 0; i < con_graph.size(); i++)
	{
		cout << "\toutput edge " << i << endl; 
		for(int j= 0; j < con_graph[i].size(); j++)
		{
			cout << "\t" << con_graph[i][j].next << endl;
		}
	}
}

int cal_quality(string &q, int begin)
{
	int sum = 0;
	int count = 0;
	for(; count < K && (count+begin) < q.size(); count++){
		sum += int(q[begin+count]);
	}
	if(count == 0)
		return 100;
	return sum/count;
}
int cal_quality_and_min(string &q, int begin, int &min)
{
	int sum = 0;
	int count = 0;
	min = 200;
	for(; count < K && (count+begin) < q.size(); count++){
		sum += int(q[begin+count]);
		min = min > int(q[begin+count])? int(q[begin+count]) : min;
	}
	if(count == 0)
		return 0;
	return sum/count;
}
/*********************************************
 *	Description(P1):
 *	
 *	Author : Zheng quangang
 *	Date   : 2015-04-10
 *********************************************/
static bool DONE = false;
void* read_thread(void *arg) {
	ARG *arg_p = (ARG *) arg;
	int threshold = arg_p->threshold;
	int minthres = arg_p->threshold;
	float percent = arg_p->percent;
	ifstream *left = arg_p->left;
	ifstream *right = arg_p->right;
	unsigned long (*buff)[BUFFER_SIZE];
	buff = arg_p->buffer;
	pthread_mutex_t *buff1 = arg_p->buffer1;
	pthread_mutex_t *buff2 = arg_p->buffer2;
	pthread_cond_t *cond1 = arg_p->cond1;
	pthread_cond_t *cond2 = arg_p->cond2;
	DeBruijnGraph *db = arg_p->db;

	int count = 0;
	int flag = 0;
	int lengths = 0;
	int lengths_cutoff = 0;
	unsigned long first_code;

	string read1, read2, read1RC, read2RC, line1, line2;

	while(!(*left).eof() && !(*right).eof())
	{
		getline(*left, line1);
		getline(*right, line2);

		if (!getline(*left, read1))
		{
			break;
		}

		getline(*right, read2RC);
		getline(*left, line1);
		getline(*right, line2);
		getline(*left, line1);
		getline(*right, line2);

		if (read1.size() < K || read2RC.size() < K)
			continue;

		read1RC = db->reverse_complement_replace_N(read1);
		read2 = db->reverse_complement_replace_N(read2RC);

		int min = 0;
		lengths = read1.size();
		lengths_cutoff = int((percent * lengths - K + 1));
		for(int index = 0; index < lengths_cutoff && index < READ_LENGTH_CUTOFF - K + 1; index ++ )
		{
			if(FILTER && (cal_quality_and_min(line1,index,min) < threshold || min < minthres)){
				continue;
			}
			first_code = 0;
			for(int i = index; i < index+K; i++)
			{
				first_code = (first_code << 2) | ((read1[i] & 0x06) >> 1);
			}
			buff[flag][++count] = first_code;
			// kmer_seq.push_a_kmer(first_code);
			first_code = 0;
			for(int i = lengths - index - K; i < lengths - index; i++)
			{
				first_code = (first_code << 2) | ((read1RC[i] & 0x06) >> 1);
			}
			// kmer_seq.push_a_kmer(first_code);
			buff[flag][++count] = first_code;
		}

		lengths = read2RC.size();
		lengths_cutoff = int((percent * lengths - K + 1));
		for(int index = 0; index < lengths_cutoff && index < READ_LENGTH_CUTOFF - K + 1; index ++ )
		{
			if(FILTER && (cal_quality_and_min(line2,index,min) < threshold || min < minthres))
				continue;
			first_code = 0;
			for(int i = lengths - index - K; i < lengths - index; i++)
			{
				first_code = (first_code << 2) | ((read2[i] & 0x06) >> 1);
			}
			// kmer_seq.push_a_kmer(first_code);
			buff[flag][++count] = first_code;
			first_code = 0;
			for(int i = index; i < index+K; i++)
			{
				first_code = (first_code << 2) | ((read2RC[i] & 0x06) >> 1);
			}
			// kmer_seq.push_a_kmer(first_code);
			buff[flag][++count] = first_code;
		}
		if(count > BUFFER_SIZE - 1000) {
			buff[flag][0] = count;
			count = 0;
			if(flag) {
				// pthread_cond_signal(cond2);
				pthread_mutex_lock(buff1); 
				pthread_mutex_unlock(buff2);
				while(buff[0][0] != 0){
					pthread_cond_wait(cond1, buff1);
				}
			}else{
				// pthread_cond_signal(cond1);
				pthread_mutex_lock(buff2); 
				pthread_mutex_unlock(buff1);
				while(buff[1][0] != 0){
					pthread_cond_wait(cond2, buff2);
				}
			}
			flag ^= 1;
		}
	}
	buff[flag][0] = count;
	count = 0;
	if(flag) {
		pthread_mutex_unlock(buff2);

	}else{
		pthread_mutex_unlock(buff1);
	}
	DONE = true;
	cout << "[INFO] read file thread end.\n";
	pthread_exit(NULL);
}
void* hash_thread(void *arg) {
	ARG *arg_p = (ARG *) arg;
	unsigned long (*buff)[BUFFER_SIZE] = arg_p->buffer;
	pthread_mutex_t *buff1 = arg_p->buffer1;
	pthread_mutex_t *buff2 = arg_p->buffer2;
	pthread_cond_t *cond1 = arg_p->cond1;
	pthread_cond_t *cond2 = arg_p->cond2;
	Kmer_Hash_Md5 *kmer_seq = arg_p->kmer_seq;

	int count = 1;
	int flag = 0;
	pthread_mutex_lock(buff1);
	count = buff[flag][0];
	bool isZero = false;
	while(1){
		if(count == 0){
			if(isZero == false){
				isZero = true;
			}else{
				if(DONE)
					break;
			}
		}
		for(int i=1; i<=count; i++){
			kmer_seq->push_a_kmer(buff[flag][i]);
		}
		if(flag) {
			buff[1][0] = 0;
			pthread_cond_signal(cond2);
			pthread_mutex_unlock(buff2);
			pthread_mutex_lock(buff1); 
			count = buff[0][0];
		}else{
			buff[0][0] = 0;
			pthread_cond_signal(cond1);
			pthread_mutex_unlock(buff1);
			pthread_mutex_lock(buff2); 
			count = buff[1][0];
		}
		flag ^= 1;
	}
	if(flag)
		pthread_mutex_unlock(buff2);
	else
		pthread_mutex_unlock(buff1);
	cout << "[INFO] hash thread end.\n";
	pthread_exit(NULL);
}
/*
 *	END(P1)
 *************************************************/
void DeBruijnGraph::initialize_kmer(Lib &lib)
{
	cout << "Building kmer hash table" << endl;
	kmer_size = K;
	int threshold = 0;
	int minthres = 0;
	unsigned long sum = 0;
	int count = 0;

	Kmer kmer1, kmer2, kmer1RC, kmer2RC;

	ifstream fin_left, fin_right;
	string left_read_file, right_read_file;
	string read1, read2, read1RC, read2RC, line1, line2;

	left_read_file = lib.get_left_read_file(0);
	right_read_file = lib.get_right_read_file(0);

	fin_left.open(left_read_file.c_str());
	fin_right.open(right_read_file.c_str());

	if(FILTER){
		while(!fin_left.eof() && !fin_right.eof())
		{
			getline(fin_left, line1);
			getline(fin_right, line2);

			if (!getline(fin_left, read1))
			{
				break;
			}

			getline(fin_right, read2RC);
			getline(fin_left, line1);
			getline(fin_right, line2);
			getline(fin_left, line1);
			getline(fin_right, line2);

			if (read1.size() < K || read2RC.size() < K)
				continue;

			for(int index = 0; index < read1.size() - K + 1 && index < READ_LENGTH_CUTOFF - K + 1; index ++ )
			{
				sum += cal_quality(line1,index);
				count ++;
			}
			for(int index = 0; index < read2RC.size() - K + 1 && index < READ_LENGTH_CUTOFF - K + 1; index ++ )
			{
				sum += cal_quality(line2,index);
				count ++;
			}
			if(count >= 1000000)
				break;
		}
	}

	float percent = 1.0;
	unsigned long first_code = 0;
	threshold = 0;
	minthres = 0;
	if(FILTER){
		threshold = round(((double)sum / count)*0.90 + 0.5);
		minthres = round(((double)sum / count)*0.80 + 0.5);
		percent = 0.95;
		cout << "[Test] threshold of kmer quality = " << threshold << endl;
		fin_right.seekg(ios::beg);
		fin_left.seekg(ios::beg);
	}
	count = 0;
	pthread_t read_t, hash_t;
	pthread_mutex_t buffer1 = PTHREAD_MUTEX_INITIALIZER, buffer2 = PTHREAD_MUTEX_INITIALIZER;
	pthread_cond_t cond1 = PTHREAD_COND_INITIALIZER, cond2 = PTHREAD_COND_INITIALIZER;
	unsigned long (*buffer)[BUFFER_SIZE] = new unsigned long[2][BUFFER_SIZE];
	buffer[0][0] = buffer[1][0] = 0;
	pthread_mutex_init(&buffer1, NULL);
	pthread_mutex_init(&buffer2, NULL);
	ARG read_struct, hash_struct;

	read_struct.threadi = read_t;
	read_struct.threshold = threshold;
	read_struct.minthres = minthres;
	read_struct.percent = percent;
	read_struct.left = &fin_left;
	read_struct.right = &fin_right;
	read_struct.buffer = buffer;
	read_struct.buffer1 = &buffer1;
	read_struct.buffer2 = &buffer2;
	read_struct.cond1 = &cond1;
	read_struct.cond2 = &cond2;
	read_struct.db = this;
	pthread_mutex_lock(&buffer1);

	if(pthread_create(&read_t, NULL, read_thread, (void*) (&read_struct) ) )
	{
		cout << "[ERROR] create read file thread error."<< endl;
		exit(1);
	}
	cout << "[INFO] read file thread has been created.\n";

	hash_struct.threadi = hash_t;
	hash_struct.buffer = buffer;
	hash_struct.buffer1 = &buffer1;
	hash_struct.buffer2 = &buffer2;
	hash_struct.cond1 = &cond1;
	hash_struct.cond2 = &cond2;
	hash_struct.db = this;
	hash_struct.kmer_seq = &kmer_seq;

	if(pthread_create(&hash_t, NULL, hash_thread, (void*) (&hash_struct) ) )
	{
		cout << "[ERROR] create hash thread error."<< endl;
		exit(1);
	}
	cout << "[INFO] hash thread has been created.\n";

	pthread_join(read_t, NULL);
	pthread_join(hash_t, NULL);

	pthread_mutex_destroy(&buffer1);
	pthread_mutex_destroy(&buffer2);

	delete[] buffer;

	fin_left.close();
	fin_right.close();

	//kmer_seq.initialize_kmer_array();
	set_lamada();
	kmer_seq.filter2();	// delete error kmers
	// kmer_seq.clear();

	Kmer_Short shortKmer;
	unsigned long code = 0;
	int num = 0;

// kmers.data
	// ifstream in(string("kmers.data").c_str());
	// while(in >> code){
	// 	in >> num;
	// 	kmer1 = code;
	// 	kmer_short.push_a_kmer_short(kmer1, num);
	// 	shortKmer = kmer1.get_next_node();
	// 	kmer_short.push_a_kmer_short(shortKmer);
	// }
ofstream out("kmers");
	KNode *p;
	while((p = kmer_seq.next_node())){
		out << p->kmer << " " << p->num << endl;
		kmer1 = p->kmer; 
		kmer_short.push_a_kmer_short(kmer1, p->num);
		shortKmer = kmer1.get_next_node();
		kmer_short.push_a_kmer_short(shortKmer);
	}
out.close();
	// kmer_seq.clear();
 // 	for (unsigned long i = 0; i < hash_length; i++)
	// {
	// 	Kmer_Node_Ptr p = kmer_seq[i];
	// 	while(p != NULL){
	// 		kmer1 = p->kmer;
	// 		kmer_short.push_a_kmer_short(kmer1, p->num);
	// 		shortKmer = kmer1.get_next_node();
	// 		kmer_short.push_a_kmer_short(shortKmer);
	// 		p = p->next;
	// 	}
	// }
	// kmer_short.initialize_kmer_short_array();
}

// void DeBruijnGraph::output_kmer()
// {
// 	cout << "Outputing all k-mer" << endl;
// 	for(unsigned int it = 0; it < kmer_seq.size(); it++)
// 	{
// 		cout << "	kmer: " << kmer_seq[it] << "	cov: " << kmer_seq.get_kmer_num(it) << endl;
// 	}
// }

/*
void DeBruijnGraph::initialize_graph()
{
	cout << "Building de Bruijn graph" << endl;
	//cout << "kmer_short.size = " << kmer_short.size() << endl;
	graph.resize(kmer_short.size());
	double count = 0;
	unsigned long hash_length = kmer_short.hash_table_length();
	for(unsigned int it = 0; it < hash_length; it++)
	{
		Kmer_Short_Node_Ptr p = kmer_short[it];
		while(p != NULL){
			graph[count].kmer_short = p->kmer_short;
			p->pos = count;
			p = p->next;
			count += 1;
		}
	}
	assert(count == kmer_short.size());

	Kmer_Short node;
	Kmer kmer;
	Kmer_Short next_node;

	for(unsigned int it = 0; it < count; it++)
	{
		node = graph[it].kmer_short;
		for(unsigned int i = 0; i < 4; i++)
		{
			kmer = node.get_next_kmer(i);
			if (kmer_seq.get_kmer_num(kmer) != 0)
			{
				next_node = kmer.get_next_node();
				graph[it].next[i] = kmer_short.has_kmer(next_node)->pos;
				graph[it].cov[i] = kmer_seq.get_kmer_num(kmer);
				assert(graph[it].cov[i] != 0);
			}
		}
	}
	
	cout << "\tnode number = " << graph.size() << endl;
	cout << "\tedge number = " << kmer_seq.size() << endl;
	kmer_seq.clear();
	kmer_short.clear();

}
*/
void DeBruijnGraph::initialize_graph()
{
	unsigned long edge_size = kmer_seq.size();
	kmer_seq.clear();
	cout << "Building de Bruijn graph" << endl;
	// graph.resize(kmer_short.size());
	cout << "[TEST] kmer_short.size() = " << kmer_short.size() << endl;
	kmer_short.esti_pos();
	kmer_short.transform2graph(graph);
	// double count = 0;
	// unsigned long hash_length = kmer_short.hash_table_length();
	// for(unsigned int it = 0; it < hash_length; it++)
	// {
	// 	Kmer_Short_Node_Ptr p = kmer_short[it];
	// 	while(p != NULL){
	// 		graph[count].kmer_short = p->kmer_short;
	// 		p->pos = count;
	// 		p = p->next;
	// 		count += 1;
	// 	}
	// }
	// assert(count == kmer_short.size());

	// Kmer_Short node;
	// Kmer kmer;
	// Kmer_Short next_node;
	// unsigned short *next;
	// unsigned long tmp_first_code;

	// for(unsigned int it = 0; it < count; it++)
	// {
	// 	node = graph[it].kmer_short;
	// 	next = kmer_short.get_kmer_nums(node);
	// 	for(unsigned int i = 0; i < 4; i++)
	// 	{
	// 	// 	kmer = node.get_next_kmer(i);
	// 	// 	cout << kmer  << " num = " << kmer.first_code << endl;
	// 	// 	// cout << kmer.get_next_node().first_code << endl;
	// 		if(next[i] == 0)
	// 			continue;
	// 		tmp_first_code = node.first_code << (64 - 2*K + 4);
	// 		tmp_first_code = tmp_first_code >> (64 - 2*K + 2);
	// 		next_node = (tmp_first_code | i);

	// 		graph[it].next[i] = kmer_short.has_kmer(next_node)->pos;
	// 		graph[it].cov[i] = next[i];
	// 		assert(graph[it].cov[i] != 0);
	// 	}
	// }
	
	cout << "\tnode number = " << graph.size() << endl;
	cout << "\tedge number = " << edge_size << endl;
	kmer_short.clear();

}

void DeBruijnGraph::output_graph(string tmp)
{
	ofstream fout(tmp.c_str());
	cout << "Outputing de Bruijn graph" << endl;
	for(int index = 0; index < graph.size(); index ++)
	{
		for (int i = 0; i < 4; i++)
		{
			if(graph[index].cov[i] > 0)
				fout << index << "\t" << graph[index].next[i] << "\t" << graph[index].cov[i]<<endl ;
		}
	}
	fout.close();
}

void DeBruijnGraph::output_condensed_de_Bruijn_graph(string file_name)
{
	cout << "Outputing condensed de Bruijn graph" << endl;
	ofstream out(file_name.c_str());

	for(int i = 0; i < con_graph.size(); i++)
	{
		for(int j = 0; j <  con_graph[i].size(); j++)
		{

			//assert(con_graph[i][j].edge_seq.size() - K + 1 == con_graph[i][j].cov.size());
			out << (i + 2) << "\t" << (con_graph[i][j].next + 2) << "\t" ;
			for(int k = 0; k < con_graph[i][j].edge_seq.size(); k ++)
			{
				out << (con_graph[i][j].edge_seq[k]);
			}
			out << endl;

			for(int k = 0; k < con_graph[i][j].cov.size(); k++)
			{
				out << con_graph[i][j].cov[k] << " ";
			}
			out << endl;
		}
	}
	out.close();
}

void DeBruijnGraph::output_de_Bruijn_graph(string file_name)
{
	cout << "Outputing de Bruijn graph" << endl;
	ofstream out(file_name.c_str());

	for(int index = 0; index < graph.size();index++)
	{
		for(int i = 0; i < 4; i++)
		{
			if(graph[index].cov[i] > 0)
				out << index << "\t" << graph[index].next[i] << "\t" << graph[index].cov[i] << endl;  
		}
	}

	out.close();
}

void DeBruijnGraph::output_min_cost_flow(string file_name)
{
//	cout << "Output minimum cost flow in math format" << endl;
	int count = 0;

	for(int i = 0; i < con_graph.size(); i++)
	{
		count += con_graph[i].size();
	}

	ofstream out(file_name.c_str());

	out.precision(2);

	// how many units of flow ?

	out << "p min " << (con_graph.size() + 2) << " " << (count*2 + con_graph.size()*2)<< endl;
	out << "n 1 " << count / 4 << endl;
	out << "n " << (con_graph.size() + 2) << " -" << count / 4 << endl;

	for(int i = 0; i < con_graph.size(); i++)
	{
		out << "a 1 " << (i + 2) << " 0 10000 0 " << endl;
	}

	for(int i = 0; i < con_graph.size(); i++)
	{
		for(int j = 0; j < con_graph[i].size(); j++)
		{
			out << "a " << (i+2) << " " << (con_graph[i][j].next + 2)
				<< " 0 " << get_d0(con_graph[i][j].cov) 
				<< " " << get_k1(con_graph[i][j].cov)
				<< endl;
			out << "a " << (i + 2) << " " << (con_graph[i][j].next + 2) 
				<< " 0 10000 " 
				<< get_k2(con_graph[i][j].cov)
				<< endl;
		}
	}

	for(int i = 0; i < con_graph.size(); i++)
	{
		out << "a " << (i + 2) << " " << (con_graph.size() + 2) << " 0 10000 0 " << endl;
	}

	out.close();

}

double DeBruijnGraph::get_d0(vector<int> &cov)
{
	double sum = 0;
	for(int i = 0; i < cov.size(); i++)
	{
		sum += cov[i];
	}
	return sum/(cov.size()*lamada);
}

double DeBruijnGraph::get_cost(vector<int> &cov, double di)
{
	double delta_1 = 0;
	double delta_2 = 0;
	for(int i = 0; i < cov.size(); i ++)
	{
		delta_1 += cov[i];
		delta_2 += n - cov[i];
	}
	double cost_di;
	cost_di = - log(di)*delta_1 - log(N - di)*delta_2;
	return cost_di;
}

double DeBruijnGraph::get_k1(vector<int> &cov)
{
	double d0 = get_d0(cov);
	double c1 = get_cost(cov, d0);
	double c0, k1;
	if(d0 > 4.0)
	{
		c0 = get_cost(cov, d0 - 2.0);
		k1 = (c1 - c0)/2.0;
	}else if( d0 > 2.0)
	{
		c0 = get_cost(cov, d0 - 1.0);
		k1 = (c1 - c0)/1.0;
	}else if( d0 > 0.5)
	{
		c0 = get_cost(cov, d0 - 0.4);
		k1 = (c1 - c0)/0.4;
	}else
	{
		c0 = get_cost(cov, d0*0.5);
		k1 = (c1 - c0)/(0.5*d0);
	}
	return k1;
}

double DeBruijnGraph::get_k2(vector<int> &cov)
{
	double d0 = get_d0(cov);
	double c1 = get_cost(cov, d0);
	double c2, k2;

	c2 = get_cost(cov, d0 + 2.0);
	k2 = (c2 - c1)/2.0;
	if (k2 < 0.0)
		cout << "Warning: k2  < 0.0" << endl;
	return k2;

}

void DeBruijnGraph::output_parameter(string file_name)
{
	ofstream fout(file_name.c_str());

	fout << "lambda=" << lamada << endl;
	fout << "K=" << K << endl;
	
	fout.close();
}

// zheng
void DeBruijnGraph::output_cdbg_com(string file_name)
{
	cout << "Outputing Cdbg and Component" << endl;

	long size = con_graph.size();
	vector<int> indegree(size,0);
	vector<int> outdegree(size,0);
	vector<vector<int> > copy_num;
	copy_num.resize(size);

	/*
	long count = 0;
	double sum = 0.0;
	for(int i=0; i < con_graph.size(); i++){
		for(int j=0; j < con_graph[i].size(); j++){
			vector<int> &tmp = con_graph[i][j].cov;
			for(int k=0; k < tmp.size(); k++){
				sum += tmp[k];
				count ++;
			}
		}
	}
	lamada = sum / count;
	cout << "\tlamada = " << lamada << endl;
	 */

	for(int i = 0; i < con_graph.size(); i++)
	{
		for(int j = 0; j <  con_graph[i].size(); j++)
		{
			long next = con_graph[i][j].next;
			indegree[next] ++;
			outdegree[i] ++;

			long sum = 0;
			int k = 0;
			for(; k < con_graph[i][j].cov.size(); k++)
			{
				sum += con_graph[i][j].cov[k];
			}
			int copy = int( round( (double(sum)/k) / lamada ) );
			copy_num[i].push_back(copy);
		}
	}

	ofstream out(file_name.c_str());
	ofstream fout("component_0");

	long ctgindex = 0;
	long comindex = 0;
	for(int i = 0; i < con_graph.size(); i++)
	{
		for(int j = 0; j <  con_graph[i].size(); j++)
		{
			long next = con_graph[i][j].next;
			int copy = copy_num[i][j];
			if(copy <= 0 && indegree[i] != 0 && outdegree[next] != 0)
				copy = 1;

			out << ">seq_" << ctgindex << "\t" << copy << "\n" ;
			for(int k = 0; k < con_graph[i][j].edge_seq.size(); k ++)
			{
				out << (con_graph[i][j].edge_seq[k]);
			}
			out << endl;

			//if(copy == 1 && con_graph[i][j].edge_seq.size() >= 2*K)
			if(copy == 1)
			{
				fout << ">component " << comindex << endl;
				fout << ctgindex << endl << endl;
				comindex ++;
				ctgindex ++;
			}else{
				ctgindex ++;
			}
		}
	}
	fout.close();
	out.close();
}
/*
void DeBruijnGraph::chop_kmer(){
	cout << "Building kmer hash from condensed graph" << endl;

	Kmer kmer;
	vector<char>::iterator it;
	for(int i = 0; i < con_graph.size(); i++)
	{
		for(int j = 0; j <  con_graph[i].size(); j++)
		{
			int len = con_graph[i][j].edge_seq.size();
			it = con_graph[i][j].edge_seq.begin();
			int k = 0;
			for(; k < len-K+1 && k < READ_LENGTH_CUTOFF - K + 1; k ++)
			{
				kmer = string(it+k,it+k+K);
				k_hash.push_a_kmer(kmer,i,j,len-k-K);
			}
			if(len > 2*READ_LENGTH_CUTOFF-K+1){
				k = len - READ_LENGTH_CUTOFF;
			}
			for(; k < len-K+1; k ++)
			{
				kmer = string(it+k,it+k+K);
				k_hash.push_a_kmer(kmer,i,j,len-k-K);
			}
		}
	}
}
void DeBruijnGraph::map_reads(Lib &lib){
	cout << "Mapping reads to condensed graph." << endl;
		kmer_size = K;

		Kmer kmer1, kmer2, kmer1RC, kmer2RC;

		ifstream fin_left, fin_right;
		string left_read_file, right_read_file;
		string read1, read2, read1RC, read2RC, line1, line2;

		left_read_file = lib.get_left_read_file(0);
		right_read_file = lib.get_right_read_file(0);

		fin_left.open(left_read_file.c_str());
		fin_right.open(right_read_file.c_str());

		long count = 0;
		long mapNumber = 0;
		while(!fin_left.eof() && !fin_right.eof())
		{
			cout << "\r Reading reads :" <<count++;

			getline(fin_left, line1);
			getline(fin_right, line2);

			if (!getline(fin_left, read1))
			{
				break;
			}

			getline(fin_right, read2RC);

			getline(fin_left, line1);
			getline(fin_right, line2);
			getline(fin_left, line1);
			getline(fin_right, line2);

			if (read1.size() < K || read2RC.size() < K)
				continue;

			read1RC = reverse_complement(read1);
			read2 = reverse_complement(read2RC);


			vector<int> a; // node_i, branch_i ,node_j, branch_j ....
			int begin,end,pos;
			int index = 0;
			kmer1 = read1.substr( 0, K);
			if(k_hash.get_kmer_pos(kmer1,begin,end,pos)){
				a.push_back(begin); // node
				a.push_back(end);	// branch
				index += pos+1;

				for(; index < read1.size() - K + 1 && index < READ_LENGTH_CUTOFF - K + 1;)
				{
					kmer1 = read1.substr( index, K);
					if(k_hash.get_kmer_pos(kmer1,begin,end,pos)){
						a.push_back(begin);
						a.push_back(end);
						index += pos+1;
					}else{
						break;
					}
				}
				int size = a.size();
				if(size>=4)
					p_hash.push_a_path( a, 0 );
//				vector<int>::iterator it = a.begin();
//				for(int i=0;i <=size-4; i++){
//					vector<int> tmp(it+i,it+i+4);
//					p_hash.push_a_path( tmp, 0 );
//					//p_hash.push_a_path( vector<int>(a.begin()+i,a.begin()+i+4), 1);
//				}
			}
			a.clear();


			//read1RC
			index = 0;
			kmer1 = read1RC.substr( 0, K);
			if(k_hash.get_kmer_pos(kmer1,begin,end,pos)){
				a.push_back(begin);
				a.push_back(end);
				index += pos+1;

				for(; index < read1RC.size() - K + 1 && index < READ_LENGTH_CUTOFF - K + 1;)
				{
					kmer1 = read1RC.substr( index, K);
					if(k_hash.get_kmer_pos(kmer1,begin,end,pos)){
						a.push_back(begin);
						a.push_back(end);
						index += pos+1;
					}else{
						break;
					}
				}
				int size = a.size();
				if(size>=4)
					p_hash.push_a_path( a, 0 );
//				vector<int>::iterator it = a.begin();
//				for(int i=0;i <=size-4; i++){
//					vector<int> tmp(it+i,it+i+4);
//					p_hash.push_a_path( tmp, 0 );
//					//p_hash.push_a_path( vector<int>(a.begin()+i,a.begin()+i+4), 1);
//				}
			}
			a.clear();


			//read2
			index = 0;
			kmer1 = read2.substr( 0, K);
			if(k_hash.get_kmer_pos(kmer1,begin,end,pos)){
				a.push_back(begin);
				a.push_back(end);
				index += pos+1;

				for(; index < read2.size() - K + 1 && index < READ_LENGTH_CUTOFF - K + 1;)
				{
					kmer1 = read2.substr( index, K);
					if(k_hash.get_kmer_pos(kmer1,begin,end,pos)){
						a.push_back(begin);
						a.push_back(end);
						index += pos+1;
					}else{
						break;
					}
				}
				int size = a.size();
				if(size>=4)
					p_hash.push_a_path( a, 0 );
//				vector<int>::iterator it = a.begin();
//				for(int i=0;i <=size-4; i++){
//					vector<int> tmp(it+i,it+i+4);
//					p_hash.push_a_path( tmp, 0 );
//					//p_hash.push_a_path( vector<int>(a.begin()+i,a.begin()+i+4), 1);
//				}
			}
			a.clear();


			//read2RC
			index = 0;
			kmer1 = read2RC.substr( 0, K);
			if(k_hash.get_kmer_pos(kmer1,begin,end,pos)){
				a.push_back(begin);
				a.push_back(end);
				index += pos+1;

				for(; index < read2RC.size() - K + 1 && index < READ_LENGTH_CUTOFF - K + 1;)
				{
					kmer1 = read2RC.substr( index, K);
					if(k_hash.get_kmer_pos(kmer1,begin,end,pos)){
						a.push_back(begin);
						a.push_back(end);
						index += pos+1;
					}else{
						break;
					}
				}
				int size = a.size();
				if(size>=4)
					p_hash.push_a_path( a, 0 );
//				vector<int>::iterator it = a.begin();
//				for(int i=0;i <=size-4; i++){
//					vector<int> tmp(it+i,it+i+4);
//					p_hash.push_a_path( tmp, 0 );
//					//p_hash.push_a_path( vector<int>(a.begin()+i,a.begin()+i+4), 1);
//				}
			}
			a.clear();
		}
		fin_left.close();
		fin_right.close();
		cout << endl << "Number of mapped reads = " << p_hash.size() << endl;
}
void DeBruijnGraph::find_paths(){
	unsigned int csize = con_graph.size();
//	map_index.resize(4*csize);

	vector<vector<long> >  in_degree;
	in_degree.resize(4*csize);
	vector<vector<long> > out_degree;
	out_degree.resize(4*csize);

	for(int i = 0; i < csize; i++)
	{
		for(int j = 0; j <  con_graph[i].size(); j++)
		{
			P_Node *tem = p_hash.has_path1(i,j,0);
			while(tem){
				if(tem->path[1] != j){
					tem = tem->next;
					continue;
				}
				vector<int>::iterator it = tem->path.begin();
				long pre_i = *(it++);
				long pre_j = *(it++);
				long i,j;
				for(;it!=tem->path.end();){
					i = *(it++);
					j = *(it++);

					in_degree[4*i+j].push_back(4*pre_i+pre_j);
					out_degree[4*pre_i+pre_j].push_back(4*i+j);

					pre_i = i;
					pre_j = j;
				}
				tem = tem->next;
			}
		}
	}


	vector<vector<long> > tmp_merge;
	tmp_merge.resize(4*csize);
	vector<long> flag;
	long count = 0;

	for(int i = 0; i < csize; i++)
	{
		for(int j = 0; j < con_graph[i].size(); j++)
		{
			long index = 4*i+j;
//			if( !((in_degree[index].size() != 1 || (in_degree[index].size() == 1 && out_degree[in_degree[index][0]].size() != 1)) && out_degree[index].size() == 1) )
//				continue;
			if( (in_degree[index].size() != 1 && out_degree[index].size() == 1) || ((in_degree[index].size() == 1 && out_degree[in_degree[index][0]].size() != 1) && out_degree[index].size() == 1) )
			{
				long next = out_degree[index][0];
				while(in_degree[next].size()==1){
					tmp_merge[index].push_back(next);
					if(out_degree[next].size() != 1)
						break;
					next = out_degree[next][0];
				}
			}
		}
	}
	ofstream ffout( string("preMergingCondensedcontigs.data").c_str());
	for(int i = 0; i < con_graph.size(); i++)
	{
		for(int j = 0; j <  con_graph[i].size(); j++)
		{
			ffout << (i + 2) << "\t" << (con_graph[i][j].next + 2) << "\t" ;
			for(int k = 0; k < con_graph[i][j].edge_seq.size(); k ++)
			{
				ffout << (con_graph[i][j].edge_seq[k]);
			}
			ffout << endl;

			for(int k = 0; k < con_graph[i][j].cov.size(); k++)
			{
				ffout << con_graph[i][j].cov[k] << " ";
			}
			ffout << endl;
			count ++;
		}
	}
	ffout.close();
	ffout.open( string("preMerging.data").c_str());
	for(int i = 0; i < tmp_merge.size(); i++)
	{
		if(tmp_merge[i].size() == 0)
			continue;
		ffout << "(" << i/4 << "," << i%4 << ")\t:\t";
		for(int j = 0; j <  tmp_merge[i].size(); j++)
		{
			ffout << "(" << tmp_merge[i][j]/4 << "," << tmp_merge[i][j]%4 << "),";
		}
		ffout << endl;
	}
	ffout.close();


	cout << "\tbegin merge" << endl;
//	vector<vector<long> >::iterator it = tmp_merge.begin();
	for(long i=0 ; i < tmp_merge.size() ; i++)
	{
		if( tmp_merge[i].size() == 0 )
			continue;
//		long index = it - tmp_merge.begin();
		Long_Node &begin = con_graph[i/4][i%4];
		for(int j = 0; j <  tmp_merge[i].size(); j++)
		{
//			cout << tmp_merge[i][j] << endl;
			Long_Node & nTmp = con_graph[tmp_merge[i][j]/4][tmp_merge[i][j]%4];
			begin.cov.insert(begin.cov.end(), nTmp.cov.begin(), nTmp.cov.end());

//			for(long i=0; i<nTmp.edge_seq.size(); i++)
//				cout << nTmp.edge_seq[i];

//			assert(distance(nTmp.edge_seq.begin()+K-1,nTmp.edge_seq.end()) >= 0);
			begin.edge_seq.reserve(begin.edge_seq.size()+nTmp.edge_seq.size()-K+1);
			begin.edge_seq.insert(begin.edge_seq.end(), nTmp.edge_seq.begin()+K-1, nTmp.edge_seq.end());

			begin.next = nTmp.next;
			flag.push_back(tmp_merge[i][j]);
		}
	}
	tmp_merge.clear();

	cout << "\tbegin erase" << endl;
	vector<long>::iterator ptr = flag.begin();
	while(ptr!=flag.end()){
		long index = *ptr;
		con_graph[index/4].erase(con_graph[index/4].begin()+index%4);
		ptr ++;
	}
	flag.clear();

	count = 0;
	ofstream fout( string("condensed_de_bruijn_graph_after_trimming.data").c_str());
	for(int i = 0; i < con_graph.size(); i++)
	{
		for(int j = 0; j <  con_graph[i].size(); j++)
		{
			fout << (i + 2) << "\t" << (con_graph[i][j].next + 2) << "\t" ;
			for(int k = 0; k < con_graph[i][j].edge_seq.size(); k ++)
			{
				fout << (con_graph[i][j].edge_seq[k]);
			}
			fout << endl;

			for(int k = 0; k < con_graph[i][j].cov.size(); k++)
			{
				fout << con_graph[i][j].cov[k] << " ";
			}
			fout << endl;
			count ++;
		}
	}
	fout.close();
	cout << "Outputing condensed de Bruijn graph" << endl;
	cout << "[TEST]Condensed con_graph edges = " << count << endl;

}


void DeBruijnGraph::condense_con_graph(){
	cout << "[TEST]Condense con_graph\n";
	long csize = con_graph.size();
	vector<int>  in_degree(csize,0);
	//in_degree.resize(csize);
	vector<int> out_degree(csize,0);
	//out_degree.resize(csize);

	for(long i = 0; i < csize; i++)
	{
		for(int j = 0; j <  con_graph[i].size(); j++)
		{
			long next = con_graph[i][j].next;
			if( next >= 0){
				in_degree[next]++;
				out_degree[i]++;
			}
		}
	}

	vector<int> new_pos(con_graph.size(),-1);

	list<int> initial_node;
	long index = 0;
	for(int i = 0; i < con_graph.size(); i ++)
	{
		if((in_degree[i] != 1 || out_degree[i] != 1) && !(in_degree[i] == 0 && out_degree[i] == 0) )
		{
			initial_node.push_back(i);
			new_pos[i] = index;
			index ++;
		}
	}

	vector<vector<Long_Node> > tmp_graph;
	tmp_graph.resize(index);

	list<int>::iterator it;
	int next = -1;
	long sum = 0;
	long i = 0;
	index = 0;
	for(it = initial_node.begin(); it != initial_node.end(); it++)
	{
		i = *it;
		for(int j = 0; j <  con_graph[i].size(); j++)
		{
			tmp_graph[index].push_back(con_graph[i][j]);
			next = con_graph[i][j].next;
			while(next >=0 &&in_degree[next] == 1 && out_degree[next] == 1)
			{
				Long_Node &tmp = tmp_graph[index][j];
				tmp.cov.insert(tmp.cov.end(),con_graph[next][0].cov.begin(),con_graph[next][0].cov.end());
				if( *(tmp.edge_seq.end()-1) != *(con_graph[next][0].edge_seq.begin()+K-2) )
					exit(1);
				//cout << "\n" << index << " " << string(tmp.edge_seq.begin(),tmp.edge_seq.end()) << endl;
				tmp.edge_seq.insert(tmp.edge_seq.end(),con_graph[next][0].edge_seq.begin()+K-1,con_graph[next][0].edge_seq.end());
				//tmp.next = con_graph[next][0].next;
				next = con_graph[next][0].next;
			}
			tmp_graph[index][j].next = new_pos[next];
			sum ++;
		}
		index ++;
	}
	sum = 0;
	ofstream fout( string("condensed_de_bruijn_graph_after_trimming.data").c_str());
	for(int i = 0; i < tmp_graph.size(); i++)
	{
		for(int j = 0; j <  tmp_graph[i].size(); j++)
		{
			//cout << "\r Output tmp_graph ";
			fout << (i + 2) << "\t" << (tmp_graph[i][j].next + 2) << "\t" ;
			for(int k = 0; k < tmp_graph[i][j].edge_seq.size(); k ++)
			{
				fout << (tmp_graph[i][j].edge_seq[k]);
			}
			fout << endl;

			for(int k = 0; k < tmp_graph[i][j].cov.size(); k++)
			{
				fout << tmp_graph[i][j].cov[k] << " ";
			}
			fout << endl;
			sum++;
		}
	}
//	for(int i = 0; i < long_tips.size(); i++)
//	{
//		fout << (i + 2) << "\t" << (long_tips[i].next + 2) << "\t" ;
//		for(int k = 0; k < long_tips[i].edge_seq.size(); k ++)
//		{
//			fout << (long_tips[i].edge_seq[k]);
//		}
//		fout << endl;
//
//		for(int k = 0; k < long_tips[i].cov.size(); k++)
//		{
//			fout << long_tips[i].cov[k] << " ";
//		}
//		fout << endl;
//		sum ++;
//	}
	fout.close();
	cout << "Outputing condensed de Bruijn graph" << endl;

	con_graph.swap(tmp_graph);
	tmp_graph.clear();
	cout << "[TEST]Condensed con_graph edges = " << sum << endl;
}

*/
