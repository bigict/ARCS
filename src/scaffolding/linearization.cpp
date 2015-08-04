#include "linearization.h"
#include "Kmer_Hash.h"
#include "Pos_Recorder.h"

#include <iostream>
#include <fstream>

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <pthread.h>
#include <sstream>

#define NI 1 //number of Insert_size
#define ND 3 //number of DELTA

extern int K ;
extern int INSERT_SIZE;
extern double DELTA;
extern int GENOME_LEN;
extern int CPU_NUM;
extern int PAIR_KMER_NUM;
extern int ITER;

static pthread_mutex_t mutex;
static pthread_mutex_t print_mutex;

extern int EDGE_CUTOFF;
extern int PAIR_KMER_CUTOFF;
extern int PAIR_READS_CUTOFF;

Pos_Recorder p_r;
vector<string> read_set1, read_set2;
Kmer_Hash pos;
Kmer_Hash pos1;
vector<unsigned int> len;

string rev_tem("TGAC");
string reverse_complement(string fw)
{
	string rc;
	rc.resize(fw.size());
	int length = fw.size();
	for (int i = length - 1; i >= 0; i--)
	{
		rc[length - 1 - i] = rev_tem[(fw.at(i) & 0x06) >> 1];
	}
	return rc;
}

void Linearization::initialize_edge(const std::string& contig_file)
{
	cout << "Reading contigs " << endl;
	ifstream in(contig_file.c_str());
	if (!in)
	{
		cout << contig_file << " open failed" << endl;
		exit(1);
	}
	string tmp_seq;
	int tmp_copy_num;

	while(!in.eof())
	{
		Edge ed;
		in >> tmp_seq;
		in >> tmp_copy_num;
		if(in.eof())
			break;
		in >> tmp_seq ;	
		ed.seq = tmp_seq;
		ed.copy_num = tmp_copy_num;
		edge.push_back(ed);
	}
	int genome_len = 0;
	for(int i = 0; i < edge.size(); i++)
	{
		genome_len += edge[i].seq.size() - K + 1;
	}
	GENOME_LEN = genome_len;
    cout << "\tcontigs number = " << edge.size() << endl;
	cout << "\ttotal contig length = " << GENOME_LEN << endl;
}

void Linearization::initialize_scaf()
{
//	cout << "begin input scaffold " << endl;
    char buff[100];
    sprintf(buff, "%s_%d", "component", ITER);
 
	ifstream in(buff);

	if (!in)
	{
		cout << buff << " open failed" << endl;
		exit(1);
	}

	string line;
	unsigned int count = 0;
	while(!in.eof())
	{
		getline(in, line);
		count ++;
	}

	component.resize(count / 3);
	gap.resize(count / 3);

	in.close();

	in.open(buff);

	char *word;
	count = 0;
	while(!in.eof())
	{
		getline(in, line);
		getline(in, line);
		word = strtok(const_cast<char*>(line.c_str()), " \t,");
		while(word)
		{
			component[count].push_back(atoi(word));
			word = strtok(NULL, " \t,");
		}
		getline(in, line);
		word = strtok(const_cast<char*>(line.c_str()), " \t,");
		while(word)
		{
			gap[count].push_back(atoi(word));
			word = strtok(NULL, " \t,");
		}
		count ++;
	}

}

void Linearization::print_edge()
{
	cout << "print edge" << endl;
	for(int i = 0; i < edge.size(); i++)
	{
		cout << "edge " << i << endl;
		cout << edge[i].seq << endl;
		cout << edge[i].copy_num << endl;
	}
	cout << "print edge end" << endl;
}

void Linearization::initialize_kmer_hash()
{
	long index = 0;
	Kmer kmer;
	int cutoff = INSERT_SIZE * 2;
	for(int i = 0; i < component.size(); i++)
	{
		if(component[i].size() == 0){
			len.push_back(0);
			continue;
		}
		index = 0;
		if(component[i].size() == 1){
			if( int(edge[component[i][0]].seq.size()) < EDGE_CUTOFF ){
				index += int(edge[component[i][0]].seq.size()) - K + 1;
				len.push_back(index + K);
				continue;
			}
			for (int k = 0; k < int(edge[component[i][0]].seq.size()) - K + 1; k ++)
			{
				if( index <= cutoff || index >= int(edge[component[i][0]].seq.size()) - K + 1 - cutoff ){
					kmer = edge[component[i][0]].seq.substr(k, K);
					pos.push_a_kmer(kmer, i, index);
				}
				index ++;
			}
		}else{
			int lens = int(edge[component[i][0]].seq.size());
			for(int j = 1; j < component[i].size(); j++)
			{
				lens += gap[i][j];
				lens += int(edge[component[i][j]].seq.size()) - K + 1;
			}

			for(int j = 0; j < component[i].size() - 1; j++)
			{
				for (int k = 0; k < int(edge[component[i][j]].seq.size()) - K + 1; k ++)
				{
					if( index <= cutoff || index >= lens - cutoff ){
						kmer = edge[component[i][j]].seq.substr(k, K);
						pos.push_a_kmer(kmer, i, index);
					}
					index ++;
				}
				index += gap[i][j];
			}
			for (int k = 0; k < int(edge[component[i][component[i].size()-1]].seq.size()) - K + 1; k ++)
			{
				if( index <= cutoff || index >= lens - cutoff ){
					kmer = edge[component[i][component[i].size()-1]].seq.substr(k, K);
					pos.push_a_kmer(kmer, i, index);
				}
				index ++;
			}
		}
		if (index == 0)
			len.push_back(edge[component[i][component[i].size()-1]].seq.size());
		else
			len.push_back(index + K);
	}
	pos.initialize_kmer_array();
	pos.free_buff();
}

void Linearization::initialize_kmer_hash_for_trainning(){
	cout << "[TEST] initialize_kmer_hash_for_trainning_insert_size" << endl;
	long index = 0;
	Kmer kmer;
	int cutoff = INSERT_SIZE * 2;
	for(int i = 0; i < component.size(); i++)
	{
		if(len[i] < cutoff){
			continue;
		}
		index = 0;
		int j = 0;
		for(; j < component[i].size() -1; j++)
		{
			for (int k = 0; k < int(edge[component[i][j]].seq.size()) - K + 1; k ++)
			{
				kmer = edge[component[i][j]].seq.substr(k, K);
				pos1.push_a_kmer(kmer, i, index);
				index ++;
			}
			index += gap[i][j];
		}
		for(; j < component[i].size() ; j++)
		{
			for (int k = 0; k < int(edge[component[i][j]].seq.size()) - K + 1; k ++)
			{
				kmer = edge[component[i][j]].seq.substr(k, K);
				pos1.push_a_kmer(kmer, i, index);
				index ++;
			}
		}
	}
	pos1.initialize_kmer_array();
	pos1.free_buff();
}

void Linearization::free_kmer_hash()
{
	pos.clear();
}

void input_read(string file_name1, string file_name2)
{
	cout << "Reading paired-reads file for mapping" << endl;

	string line1, line2, read1, read2;
	unsigned long count = 0;
	ifstream ifile1;
	ifstream ifile2;

	ifile1.open(file_name1.c_str());
	ifile2.open(file_name2.c_str());

	while (!ifile1.eof() && !ifile2.eof())
	{
		getline(ifile1, line1);
		getline(ifile2, line2);
		if (!getline(ifile1, read1))
		{
			break;
		}
		getline(ifile2, read2);
		getline(ifile1, line1);
		getline(ifile2, line2);
		getline(ifile1, line1);
		getline(ifile2, line2);

		count ++;
	}

	ifile1.close();
	ifile2.close();
	
	cout << "\tread set size = " << count << endl;

	read_set1.resize(count);
	read_set2.resize(count);

	count = 0;
	ifile1.open(file_name1.c_str());
	ifile2.open(file_name2.c_str());

	while (!ifile1.eof() && !ifile2.eof())
	{
		getline(ifile1, line1);
		getline(ifile2, line2);
		if (!getline(ifile1, read1))
		{
			break;
		}
		getline(ifile2, read2);
		getline(ifile1, line1);
		getline(ifile2, line2);
		getline(ifile1, line1);
		getline(ifile2, line2);

		read_set1[count] = read1;
		read_set2[count] = read2;
	
		count ++;
	}
	ifile1.close();
	ifile2.close();
}


void free_read()
{
	read_set1.clear();
	read_set2.clear();
}

static void *thread(void *ptr)
{
	int *a = (int*)ptr;
		int begin = (*a) * (read_set1.size() / CPU_NUM + 1);
		int end = (*a) * (read_set1.size() / CPU_NUM + 1) + read_set1.size() / CPU_NUM;
		if (end > read_set1.size())
		{
			end = read_set1.size();
		}

		pthread_mutex_lock(&print_mutex);
		cout << "thread " << *a << " begin map reads" << endl;
		pthread_mutex_unlock(&print_mutex);

		int count = 0;
		int p1_i, p1_j, p2_i, p2_j, len_tmp;
		int overlaped = 0;
		string read1, read2, read1RC, read2RC;
		Kmer kmer1, kmer2;
		pair<int, int> left, right;



		for(int index = begin; index < end; index++ )
		{
			read1 = read_set1[index];
			read2 = read_set2[index];
			if (read1.size() < K || read2.size() < K)
			{
				continue;
			}
			read2RC = reverse_complement(read2);
			read1RC = reverse_complement(read1);

			int flag = 0;
			string fstring = "";
			for(int i = 0; i < read1.size() - K + 1 && i < read2.size() - K + 1; i ++)
			{
				kmer1 = read1.substr( i, K);
				kmer2 = read2RC.substr( i, K);

				left = pos.get_kmer_pos(kmer1);
				right = pos.get_kmer_pos(kmer2);

				if(left.first != -1 && right.first != -1)
				{
					p1_i = left.first;
					p2_i = right.first;
					p1_j = left.second;
					p2_j = right.second;


					if(p1_i != p2_i )
					{
						pthread_mutex_lock(&mutex);
						PAIR_KMER_NUM ++;
						pthread_mutex_unlock(&mutex);
						stringstream t;
						string temp;
						t<<p1_i<<"|"<<p2_i<<"|"<<index<<" ";
						t>>temp;
						if(fstring.find(temp)==string::npos){
							flag = 1;
							fstring.append(temp);
						}else{
							flag = 0;
						}

						len_tmp = INSERT_SIZE + p1_j - p2_j;
						overlaped = long(len[p1_i]) - len_tmp - K + 1;
						if(overlaped > 0){
							if( overlaped > INSERT_SIZE )
								continue;
							len_tmp = long(len[p1_i]) - K + 1;
						}
						pthread_mutex_lock(&mutex);
						if(flag == 1){
							p_r.add_a_dis(p1_i, p2_i, 1, len_tmp,1);
							flag = 1;
						}
						else
							p_r.add_a_dis(p1_i, p2_i, 1, len_tmp,0);
						pthread_mutex_unlock(&mutex);




					}
				}

			}

			flag = 0;
			fstring.clear();
			for(int i = 0; i < read1.size() - K + 1 && i < read2.size() - K + 1; i++)
			{
				kmer1 = read2.substr(i, K);
				kmer2 = read1RC.substr(i, K);

				left = pos.get_kmer_pos(kmer1);
				right = pos.get_kmer_pos(kmer2);

				if(left.first != -1 && right.first != -1)
				{
					p1_i = left.first;
					p2_i = right.first;
					p1_j = left.second;
					p2_j = right.second;

				

					if(p1_i != p2_i)
					{
						pthread_mutex_lock(&mutex);
						PAIR_KMER_NUM ++;
						pthread_mutex_unlock(&mutex);
						stringstream t;
						string temp;
						t<<p1_i<<"|"<<p2_i<<"|"<<index<<" ";
						t>>temp;
						if(fstring.find(temp)==string::npos){
							flag = 1;
							fstring.append(temp);
						}else{
							flag = 0;
						}
						len_tmp = INSERT_SIZE + p1_j - p2_j;
						overlaped = long(len[p1_i]) - len_tmp - K + 1;
						if(overlaped > 0){
							if( overlaped > INSERT_SIZE )
								continue;
							len_tmp = long(len[p1_i]) - K + 1;
						}
						pthread_mutex_lock(&mutex);
						if(flag == 1){
							p_r.add_a_dis(p1_i, p2_i, 1, len_tmp,1);
							flag = 1;
						}
						else
							p_r.add_a_dis(p1_i, p2_i, 1, len_tmp,0);
						pthread_mutex_unlock(&mutex);

			
					}
				}
			}
		}
//cout << "[DONE] links_dis.data \n";
}

void Linearization::estimate_insert_size()
{
	cout << "Training mean and variance of insert size" << endl;
	//cout << "pos size " << pos.size() << endl;
	string read1, read2, read1RC, read2RC;
	Kmer kmer1, kmer2;

	//cout << "kmer size : " << K << endl;

	int count = 0;

	list<int> len_dis;
	int pi, pj;

	pair<int, int > right, left;
	//cout << "reads number " << read_set1.size() << endl;
	while (count < read_set1.size() && len_dis.size() < 1000)
	{
		read1 = read_set1[count];
		read2 = read_set2[count];

		if (read1.size() < K || read2.size() < K)
		{
			count ++;
			continue;
		}

		read2RC = reverse_complement(read2);
		read1RC = reverse_complement(read1);

		for(int i = 0; i < read1.size()-K+1 && i < read2.size() - K + 1; i++)
		{
			kmer1 = read1.substr( i, K);
			kmer2 = read2RC.substr( i, K);

			left = pos1.get_kmer_pos(kmer1);
			right = pos1.get_kmer_pos(kmer2);

			if(left.first != -1  && right.first != -1 && left.first == right.first)
			{
				if (right.second - left.second > INSERT_SIZE / 4 && right.second - left.second < 2*INSERT_SIZE) 
					len_dis.push_back(right.second - left.second);
			}

			kmer1 = read2.substr( i, K);
			kmer2 = read1RC.substr( i, K);

			left = pos1.get_kmer_pos(kmer1);
			right = pos1.get_kmer_pos(kmer2);

			if(left.first != -1 && right.first != -1 && left.first == right.first)
			{
				if (right.second - left.second > INSERT_SIZE / 4 && right.second - left.second < 2* INSERT_SIZE )
					len_dis.push_back(right.second - left.second);
			}
		}

		count ++;
	}

	ofstream out("insert_size");
	for(list<int>::iterator it = len_dis.begin(); it != len_dis.end(); it++)
	{
		out << *it << endl;
	}
	out.close();

	long sum = 0;
	int index = 0;

	for(list<int>::iterator it = len_dis.begin(); it != len_dis.end(); it++)
	{
		if ((*it) > 10)
		{
			sum += *it;
			index ++;
		}
	}

	double tmp_delta;
	int tmp_insert_size;

	//cout << "\ttraining insert size set = " << index << endl;

	if ( index > 0 )
	{
		tmp_insert_size = double(sum) / index;
	}
	else
	{
		cout << "SoftWare is not able to compute insert size." << endl;
		exit(1);
	}

	sum = 0;
	index = 0;

	for(list<int>::iterator it = len_dis.begin(); it != len_dis.end(); it++)
	{
		if ( (*it) > 10)
		{
			sum += pow(*it - tmp_insert_size, 2);
			index ++;
		}
	}

	tmp_delta = sqrt( double(sum) / index);

	INSERT_SIZE = tmp_insert_size;
	DELTA = tmp_delta;

	//p_r.set_delta(tmp_delta);
	//p_r.set_insert_size(tmp_insert_size);

	cout << "\tmean insert size = " << tmp_insert_size << endl;	
	cout << "\tvariance = " << tmp_delta << endl;

}

void Linearization::map_read(string file_name1, string file_name2)
{

	//cout << "len len" << len.size() << endl;
	input_read(file_name1, file_name2);
	initialize_kmer_hash_for_trainning();
	estimate_insert_size();
	pos1.clear();

	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&print_mutex, NULL);

	int *ptr = new int[CPU_NUM];
	pthread_t *id = new pthread_t[CPU_NUM];
	int ret;

	for (int i = 0; i < CPU_NUM; i++)
	{
		ptr[i] = i;
		ret = pthread_create(&(id[i]), NULL, thread, &(ptr[i]));
		if(ret != 0)
		{
			cout << "create pthread error\n" << endl;
			exit(1);
		}
	}

	for(int i = 0; i < CPU_NUM; i++)
	{
		pthread_join(id[i], NULL);
	}

	//p_r.set_pair_kmer_num(PAIR_KMER_NUM);
	cout << "\tmap paired kmer number = " << PAIR_KMER_NUM << endl;
	//PAIR_KMER_CUTOFF = (int)((double)PAIR_KMER_NUM / GENOME_LEN);
	cout << "\tpaired kmer cutoff = " << PAIR_KMER_CUTOFF << endl;
	pthread_mutex_destroy(&mutex);
	pthread_mutex_destroy(&print_mutex);
	
	free_read();
}

void Linearization::output_para()
{
	char buff[100];
	sprintf(buff, "scaffold_parameter_%d", ITER);
	ofstream fout(buff);
	fout << "K: " << K << endl;
	fout << "CPU: " << CPU_NUM << endl;
	fout << "EDGE_CLUSTER_NUM: " << component.size() << endl;
	fout << "PAIR_KMER_CUTOFF: " << PAIR_KMER_CUTOFF << endl;
	fout << "EDGE_CUTOFF: " << EDGE_CUTOFF << endl;
	fout << "GENOME_LEN: " << GENOME_LEN << endl;
	fout << "INSERT_SIZE : " << INSERT_SIZE << endl;
	fout << "DELTA : " << DELTA << endl;

	fout.close();
}



void Linearization::linearize(string read_file1, string read_file2)
{
	initialize_kmer_hash();
	p_r.resize(component.size());
	u_e_g.resize(component.size());
	/*
	for(int j = 0; j < len.size(); j++)
	{
		u_e_g.add_a_len(len[j]);
		p_r.add_a_len(len[j]);
	}
	*/
	p_r.initialize_mutex_array(component.size());

	map_read(read_file1, read_file2);

	p_r.destroy_mutex_array();

	free_kmer_hash();
	p_r.compute_dis(u_e_g);

	//u_e_g.input_real_dis();

	//u_e_g.initialize_gaussian();
	u_e_g.set_edge_score();
	u_e_g.output_graph("contig_arc_graph_before_remove_ambigous_arcs");
	u_e_g.output_detailed_graph("detailed_contig_arc_graph_before_remove_ambigous_arcs");
	
	u_e_g.set_prob_cutoff();
	u_e_g.remove_amb_edges();

	/* Zheng quangang
	 * date : 2015-03-30
	 */
	u_e_g.build_long_contig_graph();
	u_e_g.extend_lcg();
	/*
	 * End
	 */

	u_e_g.output_graph("contig_arc_graph_after_remove_ambigous_arcs");
	u_e_g.output_detailed_graph("detailed_contig_arc_graph_after_remove_ambigous_arcs");

	output_para();
	
	u_e_g.output_edge_len();
	u_e_g.output_lp();

}


