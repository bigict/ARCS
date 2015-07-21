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

#include <iostream>
#include "DeBruijnGraph.h"
#include <list>
#include <fstream>
#include <assert.h>
#include <cmath>
 #include <cstdlib>

extern int NUM_THREAD;
extern int K;
extern bool FILTER;
extern int READ_LENGTH_CUTOFF;
extern string rev_tem("TGAC");

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

void DeBruijnGraph::build_condensed_de_bruijn_graph()
{
	cout << "Building condensed de bruijn graph " << endl;

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

	//cout << "graph.size() = " << graph.size() << endl;

	//cout << "\tcalculate in_degree out_degree" << endl;
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

	//cout << "\tproduce initial node" << endl;
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
	unsigned int sum = 0;
	string tmp_kmer_seq;

	bool flag = false;

	for(it = initial_node.begin(); it != initial_node.end(); it++)
	{
	//	cout << "node " << i << " it " <<(*it) << endl;
		for(int j = 0; j < 4; j++)
		{

			index = *it;
			if (graph[index].cov[j] > 0)
			{
				tmp_kmer_seq = graph[index].kmer_short.get_seq();

				//assert(tmp_kmer_seq.size() == K - 1);
				
				for(int kmer_seq_i = 0; kmer_seq_i < tmp_kmer_seq.size(); kmer_seq_i ++)
				{
					tmp_str.push_back( tmp_kmer_seq[kmer_seq_i]);
				}
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
}

void DeBruijnGraph::split_identical_edge()
{
	//cout << "deal identical edge" << endl;
	vector<vector<int> > vec_pos;
	//vector<int> vec_next;
	int cur;
	int k;
	
	
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
				}
			}
		}

		vec_pos.clear();
	}
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
//	con_graph[i].erase(con_graph[i].begin() + j);
	size = con_graph.size();

	con_graph[i][j].edge_seq.swap(fir_str);
	con_graph[i][j].cov.swap(fir_cov);
	con_graph[i][j].next = size;
//	add_con_edge(i, fir_str, fir_cov, size);
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
		cout << " delete failed ... size not ok" << endl;
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
	cout << "compute all kmer number" << endl;
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
	
	for(int i = 0; i < graph.size(); i++)
	{
		in_degree[i] = 0;
		out_degree[i] = 0;				
	}
	//cout << "graph.size() = " << graph.size() << endl;
	//cout << "\tcalculate in_degree out_degree" << endl;
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
	int *new_pos = new int[graph.size()];
	for(int i = 0; i < graph.size(); i++)
	{
		new_pos[i] = -1;
	}

	//cout << "\tproduce initial node" << endl;
	int index = 0;
	for(int i = 0; i < graph.size(); i ++)
	{
		if (in_degree[i] != 1 || out_degree[i] != 1)
		{
			initial_node.push_back(i);
			new_pos[i] = index;
			index ++;
		}
	}

	//cout << "node num  =  " << index << endl;

	list<int> node_index;
	list<int> node_next_index;
	list<int> node_cov;
	unsigned int sum = 0;

	//int i = 0;
	list<int>::iterator it;
	unsigned int sum1 = 0;
	for(it = initial_node.begin(); it != initial_node.end(); it++)
	{
		for(int j = 0; j < 4; j++)
		{
			index = *it;
			if (graph[index].cov[j] > 0)
			{
		//		cout << "check\t" << index << "\t" << graph[index].next[j] << endl;
			
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
				if(in_degree[node_index.front()] == 0 || out_degree[graph[node_index.back()].next[node_next_index.back()]] == 0 )
				{
//					if(avg(node_cov) < lamada && node_index.size()+K-1 < 2*K)
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
//					delete_con_edge(node_index, node_next_index);
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
			//cout << it->edge_seq << endl;
			cout << "\t" << con_graph[i][j].next << endl;
		}
	}
}
/*********************************************
 *	Description(P1):
 *	
 *	Author : Zheng quangang
 *	Date   : 2015-04-10
 *********************************************/

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
static bool DONE = false;
void* read_thread(void *arg) {
	ARG *arg_p = (ARG *) arg;
	int threshold = arg_p->threshold;
	int minthres = arg_p->threshold;
	float percent = arg_p->percent;
	ifstream *left = arg_p->left;
	ifstream *right = arg_p->right;
	unsigned long (*buff)[BUFFER_SIZE][2];
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
	unsigned long first_code,second_code;

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
		for(int index = 0; index < int((percent * read1.size() - K + 1)) && index < READ_LENGTH_CUTOFF - K + 1; index ++ )
		{
			if(FILTER && (cal_quality_and_min(line1,index,min) < threshold || min < minthres)){
				continue;
			}
			first_code = 0;
			second_code = 0;
			for(int i = index; i < index+32; i++)
			{
				// if(read1[i] == 'N')
				// 	first_code = (first_code << 2) | (rand() & 0x3);
				// else
					first_code = (first_code << 2) | ((read1[i] & 0x06) >> 1);
			}
			for(int i = index+32; i < index+K; i++)
			{
				// if(read1[i] == 'N')
				// 	second_code = (second_code << 2) | (rand() & 0x3);
				// else
					second_code = (second_code << 2) | ((read1[i] & 0x06) >> 1);
			}
			buff[flag][++count][0] = first_code;
			buff[flag][count][1] = second_code;
			first_code = 0;
			second_code = 0;			
			for(int i = lengths - index - K; i < lengths - index - K + 32; i++)
			// for(int i = index; i < index+32; i++)
			{
				// if(read1RC[i] == 'N')
				// 	first_code = (first_code << 2) | (rand() & 0x3);
				// else
					first_code = (first_code << 2) | ((read1RC[i] & 0x06) >> 1);
			}
			for(int i = lengths - index - K + 32; i < lengths - index; i++)
			{
				// if(read1RC[i] == 'N')
				// 	second_code = (second_code << 2) | (rand() & 0x3);
				// else
					second_code = (second_code << 2) | ((read1RC[i] & 0x06) >> 1);
			}
			buff[flag][++count][0] = first_code;
			buff[flag][count][1] = second_code;
		}

		lengths = read2RC.size();
		lengths_cutoff = int((percent * lengths - K + 1));
		for(int index = 0; index < int((percent * read2.size() - K + 1)) && index < READ_LENGTH_CUTOFF - K + 1; index ++ )
		{
			if(FILTER && (cal_quality_and_min(line2,index,min) < threshold || min < minthres))
				continue;
			first_code = 0;
			second_code = 0;
			for(int i = lengths - index - K; i < lengths - index - K + 32; i++)
			{
				// if(read2[i] == 'N')
				// 	first_code = (first_code << 2) | (rand() & 0x3);
				// else
					first_code = (first_code << 2) | ((read2[i] & 0x06) >> 1);
			}
			for(int i = lengths - index - K + 32; i < lengths - index; i++)
			{
				// if(read2[i] == 'N')
				// 	second_code = (second_code << 2) | (rand() & 0x3);
				// else
					second_code = (second_code << 2) | ((read2[i] & 0x06) >> 1);
			}
			buff[flag][++count][0] = first_code;
			buff[flag][count][1] = second_code;
			first_code = 0;
			second_code = 0;
			
			for(int i = index; i < index+32; i++)
			{
				// if(read2RC[i] == 'N')
				// 	first_code = (first_code << 2) | (rand() & 0x3);
				// else
					first_code = (first_code << 2) | ((read2RC[i] & 0x06) >> 1);
			}
			
			for(int i = index+32; i < index+K; i++)
			{
				// if(read2RC[i] == 'N')
				// 	second_code = (second_code << 2) | (rand() & 0x3);
				// else
					second_code = (second_code << 2) | ((read2RC[i] & 0x06) >> 1);
			}
			buff[flag][++count][0] = first_code;
			buff[flag][count][1] = second_code;
		}
		if(count > BUFFER_SIZE - 1000) {
			buff[flag][0][0] = count;
			count = 0;
			if(flag) {
				// pthread_cond_signal(cond2);
				pthread_mutex_lock(buff1); 
				pthread_mutex_unlock(buff2);
				while(buff[0][0][0] != 0){
					// cout << "1 read wait.\n ";
					pthread_cond_wait(cond1, buff1);
				}
			}else{
				// pthread_cond_signal(cond1);
				pthread_mutex_lock(buff2); 
				pthread_mutex_unlock(buff1);
				while(buff[1][0][0] != 0){
					// cout << "0 read wait.\n ";
					pthread_cond_wait(cond2, buff2);
				}
			}
			flag ^= 1;
			// cout << flag << endl;
		}
	}
	buff[flag][0][0] = count;
	count = 0;
	if(flag) {
		// pthread_cond_signal(cond2);
		// pthread_mutex_lock(buff1); 
		pthread_mutex_unlock(buff2);
		// while(buff[0][0] != 0)
		// 	pthread_cond_wait(cond1, buff1);
		// buff[0][0] = 0;
		// pthread_cond_signal(cond1);
		// pthread_mutex_unlock(buff1);

	}else{
		pthread_mutex_unlock(buff1);
		// pthread_mutex_lock(buff2); 
		// while(buff[1][0] != 0)
		// 	pthread_cond_wait(cond2, buff2);
		// buff[1][0] = 0;
		// // pthread_cond_signal(cond2);
		// pthread_mutex_unlock(buff2);
	}
	DONE = true;
	cout << "[INFO] read file thread end.\n";
	pthread_exit(NULL);
}
void* hash_thread(void *arg) {
	ARG *arg_p = (ARG *) arg;
	unsigned long (*buff)[BUFFER_SIZE][2] = arg_p->buffer;
	pthread_mutex_t *buff1 = arg_p->buffer1;
	pthread_mutex_t *buff2 = arg_p->buffer2;
	pthread_cond_t *cond1 = arg_p->cond1;
	pthread_cond_t *cond2 = arg_p->cond2;
	Kmer_Hash_Md5 *kmer_seq = arg_p->kmer_seq;

	int count = 1;
	int flag = 0;
	pthread_mutex_lock(buff1);
	count = buff[flag][0][0];
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
			kmer_seq->push_a_kmer(buff[flag][i][0], buff[flag][i][1]);
		}
		if(flag) {
			buff[1][0][0] = 0;
			pthread_cond_signal(cond2);
			pthread_mutex_unlock(buff2);
			pthread_mutex_lock(buff1); 
			count = buff[0][0][0];
		}else{
			buff[0][0][0] = 0;
			pthread_cond_signal(cond1);
			pthread_mutex_unlock(buff1);
			pthread_mutex_lock(buff2); 
			count = buff[1][0][0];
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
	//Kmer_Short node1, node2, node1RC, node2RC;
	Kmer_Short node;

	ifstream fin_left, fin_right;
	string left_read_file, right_read_file;
	string read1, read2, read1RC, read2RC, line1, line2;

	
	left_read_file = lib.get_left_read_file(0);
	right_read_file = lib.get_right_read_file(0);
	
	cout << left_read_file << endl;
	cout << right_read_file << endl;

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
	unsigned long (*buffer)[BUFFER_SIZE][2] = new unsigned long[2][BUFFER_SIZE][2];
	buffer[0][0][0] = buffer[1][0][0] = 0;
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


	set_lamada();
	kmer_seq.filter2();
	Kmer_Short shortKmer;
	unsigned long code = 0;
	int num = 0;

	KNode *p;
	while((p = kmer_seq.next_node())){
		kmer1 = p->kmer;
		// cout << kmer1 << endl;
		kmer_short.push_a_kmer_short(kmer1, p->num);
		shortKmer = kmer1.get_next_node();
		kmer_short.push_a_kmer_short(shortKmer);
		// cout << "here" << endl;
	}

	// kmer_seq.filter();

	// unsigned long hash_length = kmer_seq.hash_table_length();
	// Kmer_Short shortKmer;
	// for (unsigned long i = 0; i < hash_length; i++)
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

	// for (unsigned long i = 0; i < kmer_seq.size(); i++)
	// {
	// 	kmer1 = kmer_seq[i];
	// 	kmer_short.push_a_kmer_short(kmer1.get_next_node());
	// 	kmer_short.push_a_kmer_short(kmer1.get_pre_node());
	// }

	// kmer_short.initialize_kmer_short_array();

}

// void DeBruijnGraph::output_kmer()
// {
// 	cout << "Output kmer : " << endl;
// 	for(unsigned int it = 0; it < kmer_seq.size(); it++)
// 	{
// 		cout << "	kmer: " << kmer_seq[it] << "	cov: " << kmer_seq.get_kmer_num(it) << endl;
// 	}
// }

void DeBruijnGraph::initialize_graph()
{
	unsigned long edge_size = kmer_seq.size();
	kmer_seq.clear();
	cout << "Building de Bruijn graph" << endl;
	kmer_short.esti_pos();
	kmer_short.transform2graph(graph);
	// graph.resize(kmer_short.size());
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

	// Kmer_Short *node;
	// Kmer kmer;
	// Kmer_Short next_node;
	// unsigned int *next;
	// unsigned long tmp_first_code;
	// unsigned long tmp_second_code;

	// for(unsigned int it = 0; it < count; it++)
	// {
	// 	node = & graph[it].kmer_short;
	// 	next = kmer_short.get_kmer_nums(*node);
	// 	for(unsigned int i = 0; i < 4; i++)
	// 	{
	// 		if(next[i] == 0)
	// 			continue;
	// 		kmer = node->get_next_kmer(i);
	// 		next_node = kmer.get_next_node();

	// 		graph[it].next[i] = kmer_short.has_kmer(next_node)->pos;
	// 		graph[it].cov[i] = next[i];
	// 		assert(graph[it].cov[i] != 0);
	// 	}
	// }
	
	cout << "\tnode number = " << graph.size() << endl;
	cout << "\tedge number = " << edge_size << endl;

	// kmer_seq.clear();
	kmer_short.clear();
}

void DeBruijnGraph::output_graph(string tmp)
{
	ofstream fout(tmp.c_str());
	cout << "Outputing de Bruijn graph" << endl;
	for(int index = 0; index < graph.size(); index ++)
	{
		//fout << graph[index].kmer_short << endl;
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
	cout << "Output condensed de Bruijn graph" << endl;
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
	//cout << "output minimum cost flow in math format" << endl;
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
