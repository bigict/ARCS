#include "uniq_edge_graph.h"

#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <limits.h>
#include <assert.h>
#include <stack>
#include <queue>
#include <stdio.h>
#include <pthread.h>

extern int K;
extern int GENOME_LEN;
extern int ITER;
extern int PAIR_KMER_CUTOFF;
extern int EDGE_CUTOFF;
extern double LINK_QUALITY_PERCENT;
extern double DELTA;
extern int CPU_NUM;
extern int INSERT_SIZE;
extern int PAIR_KMER_NUM;

extern vector<unsigned int> len;
static vector< list<Dis_Node> > con;
static vector< list<Dis_Node> > conl;
static vector<double> Gaussian;
static pthread_mutex_t print_mutex;

double log_poisson(double lam, int i)
{
	if( abs(lam - 0) <= 0.000001)
	{
		return -200000;
	}else
	{	
		double sum = 0;
		for (int index = 1; index <= i; index ++)
		{
			sum += log(lam) - log( index);
		}
		return sum -lam;
	}
}

void Uniq_Edge_Graph::initialize_gaussian()
{
	vector<double> gaussian;

	gaussian.resize(2*INSERT_SIZE);
	Gaussian.resize(2*INSERT_SIZE);
	
	double temp0, temp1, temp2;

	for(int i = 0; i < 2*INSERT_SIZE; i++)
	{

		temp0 = -pow(i - INSERT_SIZE, 2)/(2.0*DELTA*DELTA);
		temp1 = exp(temp0);
		temp2 = sqrt(2*3.141592653)*DELTA;
		gaussian[i] = temp1 / temp2;
	}

	double sum = 0;
	
	for(int i = 0; i < 2*INSERT_SIZE; i++)
	{
		sum += gaussian[i];
		Gaussian[i] = sum;
	}	
	//cout << "Gaussian sum = " << sum << endl;
}

double get_Gaussian(int i)
{
	if (i < 0 )
		return 0;
	else if ( i >= 2*INSERT_SIZE)
		return 1;
	else
		return Gaussian[i];
}

double compute_lam(int a, int b, int d)
{
	double sum = 0;

	for (int i = 0; i < len[a] - K + 1; i++)
	{	
		sum += get_Gaussian(len[a] - K + d + len[b] - K - i + 1) - get_Gaussian(len[a] - K + d - i);
	}
	return PAIR_KMER_NUM * sum / GENOME_LEN;
}

double compute_score(int a, int b, int d, int c)
{
	double sum = 0;

	for (int i = 0; i < len[a] - K + 1; i++)
	{	
		sum += get_Gaussian(len[a] - K + d + len[b] - K - i + 1) - get_Gaussian(len[a] - K + d - i);
	}
	double lam = PAIR_KMER_NUM * sum / GENOME_LEN;
	return log_poisson(lam, c);
}

static void *thread( void *ptr )
{

	int *a = (int*)ptr;
	
	pthread_mutex_lock(&print_mutex);
	cout << "thread " << *a << " begin compute edge score" << endl;
	pthread_mutex_unlock(&print_mutex);
	unsigned int size = 0;
	if(con.size() % CPU_NUM == 0)
		size = con.size() / CPU_NUM;
	else
		size = con.size() / CPU_NUM + 1;
	unsigned int begin = (*a) * size;
	unsigned int end = (*a + 1) * size;
	
	if (end > con.size())
	{
		end = con.size();
	}
	list<Dis_Node>::iterator it;
	
	for (unsigned int i = begin; i < end; i++)
	{
		for(it = con[i].begin(); it != con[i].end(); it++)
		{
			it->score = compute_score(i, it->id, it->dis - len[i] + K, it->c);
		}
	}
//	cout << "thread " << *a << " end" << endl;
}

void Uniq_Edge_Graph::set_edge_score()
{
	cout << "compute edge score" << endl;
	
	initialize_gaussian();
	
	int *ptr = new int[CPU_NUM];
	pthread_t *id = new pthread_t[CPU_NUM];
	int ret;
	
	pthread_mutex_init(&print_mutex, NULL);
	for(int i = 0; i < CPU_NUM; i++)
	{
		ptr[i] = i;
		ret = pthread_create(&(id[i]), NULL, thread, &(ptr[i]));
		if (ret != 0)
		{
			cout << "create pthread error" << endl;
			exit(1);
		}
	}
	for(int i = 0; i < CPU_NUM; i++)
	{
		pthread_join(id[i], NULL);
	}
	pthread_mutex_destroy(&print_mutex);
}


void Uniq_Edge_Graph::add_a_dis(int i, int j, int c_dis, int cc)
{
	if(i >= con.size() | j >= con.size())
	{
		cout << "i or j out of range" << endl;
		exit(0);
	}
	
	list<Dis_Node>::iterator it;
	bool find = false;
	
	for(it = con[i].begin(); it != con[i].end(); it++)
	{
		if(it->id == j)
		{
			find = true;
			it->dis = c_dis;
			it->c = cc;
			break;
		}	
	}
	if(find == false)
	{
		Dis_Node nd;
		nd.id = j;
		nd.dis = c_dis;
		nd.c = cc;
		con[i].push_back(nd);
	}
}

/*
void Uniq_Edge_Graph::add_a_dis(int i, int j, int c_dis)
{
	if(i >= con.size() | j >= con.size())
	{
		cout << "i or j out of range" << endl;
		exit(0);
	}

	list<Dis_Node>::iterator it;
	bool find = false;
	
	for(it = con[i].begin(); it != con[i].end(); it++)
	{
		if(it->id == j)
		{
			find = true;
			it->dis = c_dis;
			break;
		}	
	}
	if(find == false)
	{
		Dis_Node nd;
		nd.id = j;
		nd.dis = c_dis;
		con[i].push_back(nd);
	}
}
*/

bool Uniq_Edge_Graph::is_a_dis(int i, int j)
{
	if(i >= con.size() | j >= con.size())
	{
		return false;
	}

	list<Dis_Node>::iterator it;
	
	for(it = con[i].begin(); it != con[i].end(); it++)
	{
		if(it->id == j)
		{
			return true;
		}	
	}
	return false;
}

int Uniq_Edge_Graph::get_dis(int i, int j)
{
	if(i >= con.size() | j >= con.size())
	{
		return 0;
	}
	
	list<Dis_Node>::iterator it;
	
	for(it = con[i].begin(); it != con[i].end(); it++)
	{
		if(it->id == j)
		{
			return it->dis;
		}
	}
	return -1;
}

void Uniq_Edge_Graph::set_prob_cutoff()
{
	cout << "get probability cutoff" << endl;
	
	list<Dis_Node>::iterator it;
	unsigned int count = 0;
	
	for(int i = 0; i < con.size(); i ++)
	{
		for(it = con[i].begin();it != con[i].end(); it++)
		{
			if (it->dis >= 0)
			{
				count ++;
			}
		}
	}

	vector<double> prob_vector;
	
	prob_vector.resize(count);

	count = 0;
	ofstream fout("all_link_score");
	for(int i = 0; i < con.size(); i ++)
	{
		for(it = con[i].begin();it != con[i].end(); it++)
		{
			if (it->dis >= 0)
			{
				prob_vector[count++] = it->score;
				fout << i << "\t" << it->id << "\t" << prob_vector[count-1] << endl;
			}
		}
	}
	fout.close();

	sort(prob_vector.begin(), prob_vector.end());


	double precise_pos = count;
	int array_pos = 0;

	precise_pos = precise_pos * LINK_QUALITY_PERCENT;
	array_pos = (int)(precise_pos + 0.5);
	log_pro_threshold = prob_vector[array_pos];
	
	if (5 * con.size() < prob_vector.size())
	{
		log_pro_threshold = prob_vector[prob_vector.size() - 3 * con.size()];
	}
	cout << "log probability threshold = " << log_pro_threshold << endl;
}


void Uniq_Edge_Graph::remove_amb_edges()
{
	cout << "begin remove ambiguous links lower than " << log_pro_threshold << endl;
	unsigned int sum = 0;
	
	for(int i = 0; i < con.size(); i++)
	{
		list<Dis_Node>::iterator it = con[i].begin(); 
		while(it != con[i].end())
		{
			if (it->score <  log_pro_threshold )
			{
				it = con[i].erase(it);
				sum ++;
			}
			else
			{
				it ++;
			}
		}
	}
	cout << "\tremoved links number = " << sum << endl;
}


int Uniq_Edge_Graph::get_c(int i, int j)
{
	for(list< Dis_Node >::iterator it = con[i].begin(); it != con[i].end(); it ++)
	{
		if(it->id == j)
		{
			return it->c;
		}
	}
	return 0;
}

void Uniq_Edge_Graph::del_a_dis(int i, int j)
{
	if(i >= con.size() | j >= con.size())
	{
		cout << "out of range" << endl;
	}
	list<Dis_Node>::iterator it = con[i].begin();
	while( it != con[i].end())
	{
		if(it->id == j)
		{
			it = con[i].erase(it);
		}else
		{
			it ++;
		}
	}
}

void Uniq_Edge_Graph::output_edge_len()
{
	char buff[100];
	sprintf(buff, "edge_cluster_len_%d", ITER);
	ofstream fout(buff);
	for (int i = 0; i < len.size(); i ++)
	{
		fout << len[i] << endl;
	}
	fout.close();
}

/*
void Uniq_Edge_Graph::output_lp()
{
	cout << "output linear programing mod" << endl;
	char buff[100];
	sprintf(buff, "position_lp_%d.lp", ITER);
	ofstream fout(buff);

    list<Dis_Node>::iterator it;
    
	fout << "Minimize\nz:";

	int count = 0;
	for(int i = 0; i < con.size(); i++)
	{
        for(it = con[i].begin(); it != con[i].end(); it++)
        {
		    fout << " + E_" << i << "_" << it->id;
		    count ++;
		    if (count % 5 == 0)
			    fout << endl;
        }
	}
	fout << endl << endl;

    fout << "Subject To\n";

	int index = 1;

	for(int i = 0; i < con.size(); i++)
	{
		for(it = con[i].begin(); it != con[i].end(); it++)
		{
			fout << " con" << index ++ << " : x_" << it->id << " - x_" << i << " + e_" << i<< "_" << it->id << " = " << it->dis << endl;
		}
    }

	for(int i = 0; i < con.size(); i++)
	{
		for(it = con[i].begin(); it != con[i].end(); it++)
		{
			fout << " con" << index ++ << " : E_" << i << "_" << it->id << " + e_" << i << "_" << it->id << " >= 0" << endl;
			fout << " con" << index ++ << " : E_" << i << "_" << it->id << " - e_" << i << "_" << it->id << " >= 0" << endl;
		}
	}

    fout << "Bounds\n";
	for (int i = 0; i < con.size(); i++)
	{
		fout << " x_" << i << " free\n";
	}
	for(int i = 0; i < con.size(); i++)
	{
		for(it = con[i].begin(); it != con[i].end(); it++)
		{
			fout << " e_" << i << "_" << it->id  << " free" << endl; 
			fout << " E_" << i << "_" << it->id  << " free" << endl;
		}
	}
	fout << endl << endl << "End";


	fout.close();
}
*/

void Uniq_Edge_Graph::output_lp()
{
	cout << "output linear programing mod" << endl;
	char buff[100];
	cout << "Output math file..." << endl;

	sprintf(buff, "position_lp_%d.math", ITER);
	ofstream fout(buff);

	int rn=0;
	for (int i = 0; i < con.size(); i++)
	{
		fout << "var x_" << i << ";" << endl;
	}

	list<Dis_Node>::iterator it;
	for(int i = 0; i < con.size(); i++)
	{
		for(it = con[i].begin(); it != con[i].end(); it++)
		{
			fout << "var e_" << i << "_" << it->id  << ";" << endl; 
			fout << "var E_" << i << "_" << it->id  << ";" << endl;
			rn++;
		}
	}
	// for(int i = 0; i <= rn/10000; i++)
	// 	fout << "var r_" << i << ";" << endl;

	fout << endl;
	fout << "minimize z:  ";

	int count = 0;
	for(int i = 0; i < con.size(); i++)
	{
        for(it = con[i].begin(); it != con[i].end(); it++)
        {
		    fout << " E_" << i << "_" << it->id << " + ";
		    count ++;
		    if (count % 10 == 0)
			    fout << endl;
        }
	}
	// for(int i = 0; i <= rn/10000; i++)
	// {
	// 	fout << "r_" << i << "+";
	// 	count ++;
	// 	if (count % 10 == 0)
	// 		fout << endl;
	// }
	fout << "0;" << endl << endl;


	int index = 1;

	for(int i = 0; i < con.size(); i++)
	{
		for(it = con[i].begin(); it != con[i].end(); it++)
		{
			fout << "s.t. con" << index ++ << " : x_" << it->id << " - x_" << i << " + e_" << i<< "_" << it->id << " = " << it->dis << ";" << endl;
		}
	}

	for(int i = 0; i < con.size(); i++)
	{
		for(it = con[i].begin(); it != con[i].end(); it++)
		{
			fout << "s.t. con" << index ++ << " : E_" << i << "_" << it->id << " + e_" << i << "_" << it->id << " >= 0;" << endl;
			fout << "s.t. con" << index ++ << " : E_" << i << "_" << it->id << " - e_" << i << "_" << it->id << " >= 0;" << endl;
		}
	}

	count = 0;
	// for(int i = 0; i < con.size(); i++)
	// {
	// 	for(it = con[i].begin(); it != con[i].end(); it++)
	// 	{
	// 		if(count%10000 == 0)
	// 		{
	// 			fout << "s.t. con" << index ++ << " : E_" << i << "_" << it->id << " + ";
	// 		}else{
	// 			if((count+1)%10000 == 0)
	// 				fout << " E_" << i << "_" << it->id << " <= r_" << count/10000 << ";" << endl;
	// 			else
	// 				fout << " E_" << i << "_" << it->id << " + ";
	// 		}
	// 		count ++;
	// 	}
	// }
	// if(count%10000 != 0)
	// 	fout << " 0 <= r_" << count/10000 << ";" << endl;

	fout << endl << "end;";


	fout.close();
}

int cmp( Edge_Seq_Element es1, Edge_Seq_Element es2)
{
	return es1.pos < es2.pos;
}


void Uniq_Edge_Graph::resize(int n)
{
	con.resize(n);
}

void Uniq_Edge_Graph::output_graph(string s)
{
	cout << "begin output edge graph" << endl;
	char buff[100];
	sprintf(buff, "%s_%d", s.c_str(), ITER);
	ofstream out(buff);
	list<Dis_Node>::iterator it;
	int count=0;
	for(int i = 0; i < con.size(); i ++)
	{
		for(it = con[i].begin();it != con[i].end(); it++)
		{
			if (it->dis > 0)
			{
				out << i << "\t" << it->id << "\t" << it->dis << "\t" << it->c << "\t" << it->score << endl;
				count++;
			}
		}
	}
	cout << "[Test] out_graph size = "<<count<<endl;
	out.close();
}


void Uniq_Edge_Graph::output_detailed_graph(string s)
{
	cout << "begin output detailed scaffold graph" << endl;

	char buff[100];
	sprintf(buff, "%s_%d", s.c_str(), ITER);

	ofstream out(buff);
	list<Dis_Node>::iterator it;
	for(int i = 0; i < con.size(); i ++)
	{
		for(it = con[i].begin();it != con[i].end(); it++)
		{
			if (it->dis > 0)
			{
				out << i << "|" << len[i] << "\t" << it->id << "|" <<  len[it->id] << "\t" << it->dis << "|" << it->c << "|" << it->score << endl;
			}
		}
	}
	out.close();
}

/*
 * Zheng QuanGang
 * Date : 2014.11.03
 */
void Uniq_Edge_Graph::input_real_dis()
{
	ifstream in("real_arc");
	cout << "[Test] input real distance...\n";
	if (!in)
	{
		cout << "real_arc open failed!" << endl;
		exit(1);
	}

	string line;
	char *word;
	int from,to,dis,c;
	double score;
	int count=0,all=0;

	list<Dis_Node>::iterator it;
	for(int i=0;i<con.size();i++)
		for(it=con[i].begin();it!=con[i].end();it++)
			all++;
	cout<<"[Test] "<<all<<" con size.\n";

	int temp = con.size();
	con.clear();
	con.resize(temp);
	while(getline(in,line))
	{
		sscanf(line.c_str(), "%d\t%d\t%d\t%d\t%lf", &from, &to, &dis, &c, &score);
		if(dis<5000)
			Uniq_Edge_Graph::add_a_dis(from, to, dis, c);
		//else
		//	Uniq_Edge_Graph::add_a_dis(to, from, -dis, c);

		count ++;
	}
	cout<<"[Test] "<<count<<" links are inputed.\n";
	in.close();
}

unsigned int Uniq_Edge_Graph::next_contig(unsigned int i){
	int max = INT_MIN;
	list<Dis_Node>::iterator maxit;
	list<Dis_Node>::iterator next = con[i].begin();
	for(; next != con[i].end(); next++){
		if((*next).score > max){
			max = (*next).score;
			maxit = next;
		}
//		if(len[num] <= INSERT_SIZE)
//			continue;
	}
	if(max != INT_MIN)
		conl[i].push_back((*maxit));
	return (*maxit).id;
}

void Uniq_Edge_Graph::remove_links(unsigned int i){
	int max = 0;
	list<Dis_Node>::iterator maxit;
	list<Dis_Node>::iterator next = con[i].begin();
	for(; next != con[i].end(); next++){
		max = max < (*next).score? (*next).score : max;
	}
	max = max * 2;
	next = con[i].begin();
	list<Dis_Node>::iterator end = con[i].end();
	while( next != end ){
		if(next->score < max){
			next = con[i].erase(next);
			end = con[i].end();
		}
		else{
			next ++;
		}
	}
}

void Uniq_Edge_Graph::build_long_contig_graph(){
	int size = con.size();
//	conl.resize(size);

	int max = INT_MIN;
	list<Dis_Node>::iterator maxit;
	vector<list<Dis_Node> >::iterator it = con.begin();
//	for(; it != con.end(); it++){
//		unsigned int count = it - con.begin();
//		if(len[count] <= INSERT_SIZE)
//			continue;
//		list<Dis_Node>::iterator next = (*it).begin();
//		for(; next != (*it).end(); next++){
//			unsigned int num = (*next).id;
//			if((*next).score > max){
//				max = (*next).score;
//				maxit = next;
//			}
//			if(len[num] <= INSERT_SIZE)
//				continue;
//		}
//		if(max != INT_MIN)
//			conl[count].push_back((*maxit));
//		max = INT_MIN;
//	}
	for(; it != con.end(); it++){
		unsigned int count = it - con.begin();
		int i = 0;
		if(con[count].size() <= 1)
			continue;
//		remove_links(count);
	}

//	ofstream out("initial_long_contig_graph.data");
//	int count=0;
//	list<Dis_Node>::iterator next;
//	for(int i = 0; i < conl.size(); i ++)
//	{
//		for(next = conl[i].begin(); next != conl[i].end(); next++)
//		{
//				out << i << "|" << len[i] << "\t" << next->id << "|" << len[next->id] << "\t" << next->dis << "|" << next->c << "|" << next->score << endl;
//				count++;
//		}
//	}
//	out.close();

//	vector< vector<unsigned int> > indegree;
//	indegree.resize(size);
//	it = conl.begin();
//	for(; it != conl.end(); it++){
//		unsigned int count = it - conl.begin();
//		list<Dis_Node>::iterator next = (*it).begin();
//		for(; next != (*it).end(); next++){
//			indegree[next->id].push_back(count);
//		}
//	}
//	for(int i=0; i<size; i++){
//		if(indegree[i].size() <= 1)
//			continue;
//		vector<unsigned int>::iterator pre = indegree[i].begin();
//		max = conl[(*pre)].begin()->score;
//		for(pre++ ; pre != indegree[i].end(); pre++){
//			if(conl[(*pre)].begin()->score > max)
//				max = conl[(*pre)].begin()->score;
//		}
//		pre = indegree[i].begin();
//		for( ; pre != indegree[i].end(); pre++){
//			if(conl[(*pre)].begin()->score < max)
//				conl[(*pre)].clear();
//		}
//	}
//
//	con.swap(conl);
}
void Uniq_Edge_Graph::extend_lcg(){


}
