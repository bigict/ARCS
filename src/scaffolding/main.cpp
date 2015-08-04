#include "Linearization.h"
#include "Lib.h"

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffolding.main"));

using namespace std;

extern int K;
extern int PAIR_KMER_CUTOFF;
extern int PAIR_READS_CUTOFF;
extern int EDGE_CUTOFF;
extern int PAIR_KMER_NUM ;
extern int GENOME_LEN;
extern int INSERT_SIZE ;
extern double DELTA ;
extern int CPU_NUM;
extern double LINK_QUALITY_PERCENT ;
extern int LOG_PRO_THRESHOLD ;
extern int ITER;

int K = 31;
int PAIR_KMER_CUTOFF = 10;
int PAIR_READS_CUTOFF = 5;
int EDGE_CUTOFF = 0;
int PAIR_KMER_NUM = 0;
int GENOME_LEN = 4000000;
int INSERT_SIZE = 180;
double DELTA = 5.0;
int CPU_NUM = 8;
double LINK_QUALITY_PERCENT = 0.01;
int LOG_PRO_THRESHOLD = -200;
int ITER = 0;

static const char *optString = "d:K:c:e:L:D:1:2:P:p:i:r:R:";

int main(int argc, char **argv) {

	cout<<"======================================================================================\n";
	cout<<"Step 3. Scaffolding\n";
	cout<<"======================================================================================\n\n";
	time_t start, end;
	start = time(NULL);

	string contig_file, scaffold_file, work_dir, read1, read2;
	int opt = 0;
	opt = getopt(argc, argv, optString);
	while(opt != -1)
	{
		switch(opt)
		{
			case 'c':
				contig_file = string(optarg);
				//cout << "contig file : " << contig_file << endl;
				break;
			//case 'l':
			//	scaffold_file = string(optarg);
			//	break;
			case 'K':
				K = atoi(optarg);
				if (EDGE_CUTOFF == 0)
				{
					EDGE_CUTOFF = K;
				}
				break;
			case 'e':
				EDGE_CUTOFF = atoi(optarg);
				if (EDGE_CUTOFF < 0)
				{
					cout << "edge length cutoff must be positive." << endl;
					exit(0);
				}
				break;
			case 'd':
				work_dir = string(optarg);
				break;
			case 'L':
				INSERT_SIZE = atoi(optarg);
				break;
			case '1':
				read1 = string(optarg);
				break;
			case '2':
				read2 = string(optarg);
				break;
			case 'P':
				LINK_QUALITY_PERCENT = atof(optarg);
				break;
			case 'p':
				CPU_NUM = atoi(optarg);
				break;
			case 'i':
				ITER = atoi(optarg);
				break;
			case 'r':
				PAIR_KMER_CUTOFF = atoi(optarg);
				break;
			case 'R':
				PAIR_READS_CUTOFF = atoi(optarg);
				break;
			default:
				cout << "incorrect parameter." << endl;
				break;
		}
		opt = getopt(argc, argv, optString);
	}

	if(chdir(work_dir.c_str()))
	{
		cout << "change work directory failed" << endl;
		exit(0);
	}

	Linearization li;
	li.initialize_edge(contig_file);
	li.initialize_scaf();
	li.linearize(read1, read2);
	
	end = time(NULL);
	cout << "scaffolding time = " << difftime(end, start) << " seconds" << endl;
	return 0;

}

