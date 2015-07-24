#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <fstream>

#include "Pos_Recorder.h"

extern int PAIR_KMER_CUTOFF;
extern int PAIR_READS_CUTOFF;
extern int EDGE_CUTOFF;
extern int K;
extern vector<unsigned int> len;
extern int INSERT_SIZE;

void Pos_Recorder::add_a_dis(int i, int j, int num, int c_dis, int reads)
{
	if(i >= con.size() | j >= con.size() )
	{
		exit(0);
		cout << "out of range" << endl;
	}

	pthread_mutex_lock(&(mutex_array[i]));

	list<Pos_Node>::iterator it;
	bool find = false;
	for(it = con[i].begin(); it != con[i].end(); it++)
	{
		if(it->id == j)
		{
			find = true;
			it->dis_sum += num * c_dis;
			it->c += num;
			it->reads += reads;
			break;
		}			
	}
	if(find == false)
	{
		Pos_Node nd;
		nd.id = j;
		nd.dis_sum = num * c_dis;
		nd.c = num;
		nd.reads = 1;
		con[i].push_back(nd);
	}

	pthread_mutex_unlock(&(mutex_array[i]));
}

void Pos_Recorder::compute_dis(Uniq_Edge_Graph &uniq_edge_graph)
{
	//int thr = 2;
	cout << "Pos Recorder computer dis to Uniq edge Graph" << endl;
	list<Pos_Node>::iterator it;
	int avg_dis;
	uniq_edge_graph.resize(con.size());
	//cout << "pos con size = " << con.size() << endl;
	
	unsigned int num = 0;
	ofstream fout("links_paired_reads");
	ofstream out("links_dis.data");
	for(int i = 0; i < con.size(); i++)
	{
		int j = 0;
		for(it = con[i].begin(); it != con[i].end(); it++,j++)
		{
			avg_dis = it->dis_sum / it->c;
			fout<<i<<"\t"<<it->id<<"\t"<<it->c<<"\t"<<it->reads<<endl;
			if(it->c >= PAIR_KMER_CUTOFF && it->reads >= PAIR_READS_CUTOFF)
			{
				num ++;
				uniq_edge_graph.add_a_dis(i, it->id, avg_dis, it->c);
				out << i << "\t" << j<< "\t" << avg_dis << endl;
			}
		}
	}
	cout << "[Test] Pos_Recorder::compute_dis : it->c >= " << PAIR_KMER_CUTOFF << " && it->reads >= " << PAIR_READS_CUTOFF << endl;
	cout << "total links num = " << num << endl;
out.close();
}

void Pos_Recorder::resize(int n)
{
	con.resize(n);
}

void Pos_Recorder::initialize_mutex_array(int n)
{
	mutex_array.resize(n);
	for(int i = 0; i < n; i ++)
	{
		pthread_mutex_init(&(mutex_array[i]), NULL);
	}
}

void Pos_Recorder::destroy_mutex_array()
{
	for (int i = 0; i < mutex_array.size(); i++)
	{
		pthread_mutex_destroy(&(mutex_array[i]));
	}
}


