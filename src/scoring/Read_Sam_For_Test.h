#ifndef _READS_SAM_FOR_TEST_H
#define _READS_SAM_FOR_TEST_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <numeric>
#include <algorithm>

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Para.h"
#include "Constant.h"
#include "Count_Struct.h"
#include "Split_String.h"
#include "Read_Sam_For_Training.h"

using namespace std;

struct Struct_Sam;

extern int overlap;
extern int var;
extern int MU;
extern int K;
extern string TEMP_FILE_DIR;

class Read_Sam_For_Test
{

friend ostream &operator<<(ostream &, const Read_Sam_For_Test &);

private:

	string sam_file_name;                 //the sam file of scoring contigs
	string contig_file_name;              //scoring contigs file name
	
    int gap_count;                        //scoring gap count
	vector<int> gap_candi_count;          //candidate number of each gap

	vector< vector<double> > gap_candidates_likelihood;   //likelihood of gap candidates 
	vector< vector<long> > read_mapped_count;             //mapped read count of each candidate

	vector<vector<string> > all_contigs;                  //for scoring contigs

    //temp vector
	vector<string > sam_info_1;
	vector<string > sam_info_2;

    //temp vector for sequence alignment
	vector<vector<int> > matrix;
	vector<vector<int> > back_track;

public:
	
    //ctor
	explicit Read_Sam_For_Test(const string&, const string& , const string&);
	Read_Sam_For_Test();
	~Read_Sam_For_Test();
	
	void scoring(Para &); 
	void output_final_scaffolds();

private:
	//void init();
	void read_contig_file();                                      //read scoring contig file
	void filter_sam_file(const string&, const string&);           //filter the sam file
	void add_missing_info(const int, const int);                  //add missing information of PerM
	void set_two_quality(const int,const int);
	void find_max_common_seq1(const int,const int);
	void find_max_common_seq2(const int,const int);
	
	double compute_error_prob1(const Para&);
	double compute_error_prob2(const Para&);
	double compute_insert_prob(const Para&);
	
	double compute_likelihood( Para &, const int, const int);
	double compute_locate_prob( const long, Para&, const int , const int);
};

#endif /*_READS_SAM_FOR_TEST_H*/
