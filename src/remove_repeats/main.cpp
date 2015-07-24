#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include "Uniq_Edge_Graph.h"
#include <time.h>

using namespace std;

extern int K ;
extern int MAX_OVERLAP ;
extern int EDGE_NUM ;
extern int ITER;

int K = 31;
int MAX_OVERLAP = 200;
int EDGE_NUM = 10000;
int ITER = 0;

static const char *optString = "K:d:O:i:";

int main(int argc, char **argv)
{
	time_t start, end;
	start = time(NULL);

	string work_dir;	
	int opt = 0;
	opt = getopt(argc, argv, optString);
	while (opt != -1)
	{
		switch (opt)
		{	
		case 'd':
			work_dir = string(optarg);
			//cout << "work directory " << work_dir << endl;
			break;
		case 'O':
			MAX_OVERLAP = atoi(optarg);
			//cout << "max overlap " << MAX_OVERLAP << endl;
			break;
		case 'K':
			K = atoi(optarg);
			//cout << "Kmer size = " << K << endl;
			break;
		case 'i':
			ITER = atoi(optarg);
			break;
		default:
			/* You won't actually get here. */
			//cout << "no \"" << opt << "\" opinion !!" << endl;
			break;
		}
		opt = getopt(argc, argv, optString);
	}

	if(chdir(work_dir.c_str()))
	{
		cout << "change work directory failed" << endl;
		exit(0);
	}

	cout << "begin to remove repeats" << endl;
	Uniq_Edge_Graph u_e_g;
		
	u_e_g.input_parameter();
	u_e_g.input_edge_len();
	u_e_g.input_edge_pos();
	u_e_g.resize(EDGE_NUM);
	
	cout << "Edge num = " << EDGE_NUM << endl;

	u_e_g.input_edge_link("contig_arc_graph_after_remove_ambigous_arcs");
	u_e_g.input_inner_component();
	u_e_g.linearize();
	
	end = time(NULL);

	cout << "removing repeats time = " << difftime(end, start) << " seconds" << endl;
	return 0;
}
