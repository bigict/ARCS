#include "Gap_Filling.h"

#include <time.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <sstream>
#include "Constant.h"

using namespace std;


extern string EDGE_FILE_NAME;
extern string SCAFFOLD_FILE_NAME;
extern string TEMP_FILE_DIR;
extern string INITIAL_EDGE_FILE_NAME;

extern int EXTEND;
extern int STEP;
extern int K;
extern int overlap;
extern int MU;
extern int var;

int EXTEND = -1;
int STEP = -1;
int K = -1;
int overlap = -1;
int MU = -1;
int var = -1;

//string EDGE_FILE_NAME;
//string SCAFFOLD_FILE_NAME;
//string TEMP_FILE_DIR = "./ARCS_TEMP/"; 
//string INITIAL_EDGE_FILE_NAME;



static const char *opt_str = "K:O:c:l:s:i:d:";

void print_help()
{
    cout << "USAGE : ./gap_filling -d work_directory -K kmer_size -O overlap_for_merge -c condensed_contig_file -i initial_condensed_contig_file -l line_component_file -s configure_file" << endl;
	cout << "gap_filling..." << endl;
	cout << "-K\t<int>\t\tkmer size" << endl;
	cout << "-d\t<string>\tworking directory" << endl;
	cout << "-O\t<int>\t\toverlap for merge two adjacent contigs" << endl;
	cout << "-c\t<string>\tcontig file name (fasta format) " << endl;
	cout << "-l\t<string>\tcomponent file name " << endl;
	cout << "-i\t<string>\tinitial contig file name" << endl;
    cout << "-s\t<string>\tconfigure file name" << endl;
	cout << endl;
}

inline bool is_number(const string &str)
{
	if(str.empty()) return false;
    bool dot=false;
	for(int i=0;i<str.size();++i)
	{
		if(!isdigit(str[i]) && str[i] != '.') return false;
        if(str[i] == '.') 
        {
            if(dot) return false;
            dot = true;
        }
	}
	return true;
}

int main(int argc, char ** argv)
{
	time_t start, end;;
	double seconds;

	time(&start);
		
	string initial_contig_file_name;
	string scaffold_file_name;
	string condensed_contig_file_name;
	string work_dir;
	string config_file_name;
	
	int opt;
	opt = getopt(argc, argv, opt_str);
	//cout << "argc : " << argc << endl;
	if (argc != 15 )
	{
		print_help();
		exit(1);
	}
	while(opt != -1)
	{
		switch(opt)
		{
			case 'K':
				K = atoi(optarg);
				break;
			case 'O':
				overlap = atoi(optarg);
				break;
			case 'd':
				work_dir = optarg;
				break;
			case 's':
				config_file_name = optarg;
				break;
			case 'c':
				condensed_contig_file_name = optarg;
				break;
			case 'l':
				scaffold_file_name = optarg;
				break;
			case 'i':
				initial_contig_file_name = optarg;
				break;
			default:
				print_help();
				exit(1);
		}
		opt = getopt(argc, argv, opt_str);
	}


	if(chdir(work_dir.c_str()))
	{
		cout << "No such file or directory!! Change work directory failed!!!" << endl;
		exit(EXIT_FAILURE);
	}

	ifstream fin(config_file_name.c_str());
	if (!fin)
	{
		cerr << "[Info] " << config_file_name << "No such file or directory!!!" << endl;
		exit(EXIT_FAILURE);
	}

	string line, key, value;

	vector<int> mus;
	vector<int> vars;

	size_t left,right,split;

	while(getline(fin,line))
	{
		//check annotation
		split = line.find_first_not_of("\t ");
		if(split==string::npos || line[split] == '#')
			continue;

		split = line.find_first_of("=:");
		if(split == string::npos)
			continue;
		
		key.assign(line.begin(), line.begin()+split);
		value.assign(line.begin()+split+1,line.end());

		left = key.find_first_not_of(" \t");
		right = key.find_last_not_of(" \t");
		key.assign(key.begin()+left, key.begin()+right+1);

		left = value.find_first_not_of(" \t");
		right = value.find_last_not_of(" \t");
		value.assign(value.begin()+left, value.begin()+right+1);

		if(key.substr(0,string("INSERT_SIZE").size()) == "INSERT_SIZE")
		{
			if(!is_number(value)) 
			{
				cerr << "Illegal parameter...\n" << line << endl;
				exit(1);
			}
			mus.push_back(atoi(value.c_str()));
		}
		if(key.substr(0,string("DELTA").size()) == "DELTA")
		{
			if(!is_number(value)) 
			{
				cerr << "[Info] Illegal parameter...\n" << line << endl;
				exit(1);
			}
			vars.push_back(atoi(value.c_str()));
		}
		
	}
	
    MU = *max_element(mus.begin(),mus.end());
    var = *max_element(vars.begin(), vars.end());

	STEP = MU + 3 * var;
	EXTEND = 200;
	
	Gap_Filling gf(condensed_contig_file_name, scaffold_file_name,initial_contig_file_name);

	gf.gap_filling();

    char buf[10];
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    string fileName;
    fileName += string(buf);
    fileName += string("mer");
    fileName += ".gap_filling_info";


	ofstream outfile(fileName.c_str());
	if (!outfile)
	{
		cerr << "[Info] Create " << fileName << "error!!!" << endl;
		exit(EXIT_FAILURE);
	}
	outfile << gf;

	time(&end);
	seconds = difftime(end,start);
	cout << "[Info] Running time = " << seconds << " seconds " << endl;
	
	exit(EXIT_SUCCESS);
}

