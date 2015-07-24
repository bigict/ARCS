#include "Read_Sam_For_Training.h"
#include "Read_Sam_For_Test.h"

#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

using namespace std;


extern int var;
extern int MU;
extern int K;

static const char *opt_str = "K:1:2:t:g:d:s:";

int var = -1;
int MU = -1;
int K = -1;

void print_help()
{							    
	cout << "Begin Scoring .... "<< endl;
    cout << "Usage : ./scoring -d workDirectory -K kmerSize -O MergeTwoEdgesOverlap -1 readFile1 -2 readFile2 -t UniqueEdgeForTraingParameter -g scoringContigs" << endl;
	cout << "./score : " << endl;
	cout << "-K\t<int>\t[kmer size]" << endl;
	cout << "-d\t<string>\t[working directory]" << endl;
	cout << "-s\t[parameter file]" << endl;
	cout << "-1\t<int>\t[paired read file1]" << endl;
	cout << "-2\t<int>\t[paired read file2]" << endl;
	cout << "-t\t<string>\t[unique contigs for training parameter file name] " << endl;
	cout << "-g\t<string>\t[multi-candidates sequences]" << endl;
	cout << endl;
}

int main(int argc, char ** argv)
{
	time_t start, end;
	double seconds;
	time(&start);

	string read_file_name1 ;
	string read_file_name2;
	string training_file_name ;
	string multi_seqs_file_name;
	string config_file_name;
	string work_dir;

	int opt;
	opt = getopt(argc, argv, opt_str);
	if(argc != 15 )
	{
		print_help();
		exit(1);
	}
	cout << argc << endl;

	while(opt != -1)
	{
		//cout << char(opt) << endl;
		switch(opt)
		{
			case 'K': 
				K = atoi(optarg); 
				break;
			case 's':
				config_file_name = optarg;
				break;
			case 'd':
				//work_dir = atoi(optarg); 
				work_dir = optarg; 
				break;
			case '1':
				read_file_name1 = optarg;
				break;
			case '2':
				read_file_name2 = optarg;
				break;
			case 't':
				training_file_name = optarg;
				break;
			case 'g':
				multi_seqs_file_name = optarg;
				break;
			default:
				print_help();
				exit(1);
		}
		opt = getopt(argc, argv, opt_str);
	}
	cout << work_dir << endl;

/*
	//call PerM and map contigs
	string command = "./ARCS_PerM.sh";
	command += " -1 ";
	command += read_file_name1;
	command += " -2 ";
	command += read_file_name2;
	command += " -t ";
	command += training_file_name;
	command += " -g ";
	command += multi_seqs_file_name;

	pid_t status;

	status = system(command.c_str());

	if (-1 == status)
	{
		cout << "system error!!" << endl;
		exit(1);
	}
	else
	{
		printf("exit status value = [0x%x]\n" , status);
		if (WIFEXITED(status))
		{
			if (0 == WEXITSTATUS(status))
			{
				printf("run shell for PerM successfully..\n");
			}
			else
			{
				printf("runt shell for PerM failed, script exit code : %d\n", status);
			}
		}
		else
		{
			printf("exit status = [%d]\n",WEXITSTATUS(status));
		}
	}
*/
		
	if(chdir(work_dir.c_str()))
	{
		cout << "No such file or directory!! Change work directory failed!!!" << endl;
		exit(EXIT_FAILURE);
	}
	
	ifstream fin(config_file_name.c_str());
	if (!fin)
	{
		cout << "configure file does not exist!!!" << endl;
		exit(EXIT_FAILURE);
	}
	cout << K << endl;

	//read configure file and get MU and var
	string line;
	stringstream ss;
	string key, value;
	while(getline(fin,line))
	{
		ss.clear();
		ss << line;
		ss >> key;
		ss >> value;
		cout << key << " " << value << endl;
		if (key == string("INSERT_SIZE:"))
		{
			MU = atoi(value.c_str());
		}
		if (key == string("DELTA:"))
		{
			var = atoi(value.c_str());
		}
	}
	
    cout << "insert size : " << MU << endl;
    cout << "insert size std var : " << var << endl;


	string training_sam_file_name = training_file_name + ".sam";

	Read_Sam_For_Training training(training_sam_file_name, training_file_name);
	
    char buf[10];
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    string fileName;
    fileName += string(buf);
    fileName += string("mer");
    fileName += ".parameter_scoring";


	Para para;
	training.update_parameter(para);

	ofstream out(fileName.c_str());
	out << para << endl;
	out.close();
	out.clear();

	time(&end);
	seconds = difftime(end, start);

	cout << "Training parameter for scoring end..." << endl;
	cout << "time : " <<  seconds << "  seconds " << endl;
	time(&start);

    //scoring contigs file
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    fileName.clear();
    fileName += string(buf);
    fileName += string("mer");
    fileName += ".multi_gap_count_to_score";

    //string gap_count_file_name = "multi_gap_count_to_score";
	string gap_count_file_name = fileName;

	//sam_file_name
	string scoring_sam_file_name = multi_seqs_file_name + ".sam";
	
	Read_Sam_For_Test Test(gap_count_file_name, scoring_sam_file_name, multi_seqs_file_name);
	
	Test.scoring(para);

    buf[0] = '\0';
    sprintf(buf, "%d",K);
    fileName.clear();
    fileName += string(buf);
    fileName += string("mer");
    fileName += ".scoring_likelihood";

	out.open(fileName.c_str());
	out << Test <<  endl;
	out.close();
	out.clear();

	time(&end);
	seconds = difftime(end, start);
	cout << "Soring time = " <<  seconds << "  seconds " << endl;

	return 0;
}
