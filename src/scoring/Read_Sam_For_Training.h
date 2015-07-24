#ifndef _READ_SAM_FOR_TRAINING_H
#define _READ_SAM_FOR_TRAINING_H

#include "Para.h"
#include "Constant.h"
#include "Count_Struct.h"
#include "Split_String.h"
#include "Func.h"

#include <algorithm>
#include <numeric>

extern int overlap;
extern int var;
extern int MU;
extern int K;
extern string TEMP_FILE_DIR;

class Read_Sam_For_Training
{

	friend ostream &operator<<(ostream &, const Read_Sam_For_Training&);

public:

	Read_Sam_For_Training(const string&,const string&);
	Read_Sam_For_Training();
	
	~Read_Sam_For_Training();

	void update_parameter(Para &); 

private:

	string contig_file_name;
	string sam_file_name;

	Count_Struct count_struct;

	//for filter
	vector<string > sam_info_1;
	vector<string > sam_info_2;
	
	//for alignment
	vector<vector<int> > matrix;
	vector<vector<int> > back_track;

private:
	
	void filter_sam_file(const string &, const string &);
	void update_insert( long );
	void update_error_type(vector<string> &);

	void set_parameter(Para &);
	
	void add_missing_info();
	void set_two_quality();

	// for count_struct
	void init_insert_counts(long _max_insert_len = MAX_INSERT_LEN);
	void init_error_types(int);
	string reverse(const string &);
	void read_contig_file();
	void find_max_common_seq1();
	void find_max_common_seq2();
};


#endif
