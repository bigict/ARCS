#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>

#include "DeBruijnGraph.h"
extern int K;
extern int READ_LENGTH_CUTOFF ;
extern bool FILTER ;
extern int NUM_THREAD;

int K = 65;
int NUM_THREAD = 1;
int READ_LENGTH_CUTOFF = 100;
bool FILTER = false;

static const char *optString = "s:d:K:E:p:";

int main(int argc, char **argv)
{
	cout<<"======================================================================================\n";
	cout<<"Step 1. Contiging\n";
	cout<<"======================================================================================\n\n";
	time_t start, end;
	start = time(NULL);
	cout << "generate contigs" << endl;
	string work_dir;
	string conf_name;
	string file_name_1;
	string file_name_2;
	
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
		case 's':
			conf_name = string(optarg);
			//cout << "configure file " << conf_name << endl;
			break;
		case 'E':
			FILTER = atoi(optarg);
			break;
		case 'p':
			NUM_THREAD = atoi(optarg);
			break;
		case 'K':
			K = atoi(optarg);
			if (K <= 64 || K > 96)
			{
				cout << "K must be in [65, 96]" << endl;
				exit(0);
			}
			//cout << "Kmer size = " << K << endl;
			break;
		default:
			/* You won't actually get here. */
			cout << "no \"" << opt << "\" opinion !!" << endl;
			break;
		}
		opt = getopt(argc, argv, optString);
	}
	
	Lib lib;

	string line;
	string key, value;
	ifstream fin(conf_name.c_str());
	if(!fin)
	{
		cout << "configure file does not exist" << endl;
		exit(1);
	}
	
	int first_i = 0;
	int eql_pos = 0;
	int left, right;
	while (getline(fin, line))
	{
		eql_pos = 0;
		for (first_i = 0; first_i < line.size(); first_i ++)
		{
			if (line[first_i] != ' ')
				break;
		}
			
		if (first_i < line.size()&& line[first_i] == '#')
		{
			continue;
		}

		for(int i = 0; i < line.size(); i++)
		{
			if (line[i] == '=' || line[i] == ':')
			{
				eql_pos = i;
				break;
			}
		}
		if( eql_pos == 0)
		{
			continue;
		}


		key = line.substr(0, eql_pos);
		value = line.substr(eql_pos+1);
		
		for(int i = 0; i < key.size(); i++)
		{
			if (key[i] != ' ')
			{
				left = i;
				break;
			}
		}

		for (int i = key.size() - 1; i >= 0; i --)
		{
			if (key[i] != ' ')
			{
				right = i;
				break;
			}
		}

		key = key.substr(left, right-left+1);

		for(int i = 0; i < value.size(); i++)
		{
			if (value[i] != ' ')
			{
				left = i;
				break;
			}
		}

		for (int i = value.size() - 1; i >= 0; i --)
		{
			if (value[i] != ' ')
			{
				right = i;
				break;
			}
		}

		value = value.substr(left, right-left+1);

		if (key == string("q1"))
		{
			file_name_1 = value;
			continue;
		}
		if (key == string("q2"))
		{
			file_name_2 = value;
			lib.push_two_read_file(file_name_1, file_name_2);
			continue;
		}
		if (key == string("READ_LENGTH_CUTOFF"))
		{
			READ_LENGTH_CUTOFF = atoi(value.c_str());
		}
	}
	
	if(chdir(work_dir.c_str()))
	{
		cout << "change work directory failed" << endl;
		exit(1);
	}

	//input reads
	
	//build de Bruijn graph
	srand((unsigned int)time(0));
	DeBruijnGraph dbg;
	dbg.initialize_kmer(lib);
	// dbg.set_lamada();//lambda refers to mean coverage of unique kmer.
	dbg.initialize_graph();//in
	dbg.set_n();//n refers to all kmer numbers.
	dbg.output_parameter("contig_parameter");
	dbg.build_condensed_de_bruijn_graph();//convert de Bruijn graph to condensed de Bruijn graph.
	dbg.output_condensed_de_Bruijn_graph("condensed_de_bruijn_graph_before_trimming.data");
	dbg.trim_bubble_tip();
	dbg.trim_bubble_tip();
	
	dbg.build_condensed_de_bruijn_graph();
	dbg.split_identical_edge();
	dbg.output_condensed_de_Bruijn_graph("condensed_de_bruijn_graph_after_trimming.data");
	dbg.output_min_cost_flow("min_cost_flow.DIMACS");

	end = time (NULL);

	cout << "contiging time = " << difftime(end, start) << " seconds " << endl;

	return 0;

}
