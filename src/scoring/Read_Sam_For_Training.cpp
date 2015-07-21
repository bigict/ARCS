#include "Read_Sam_For_Training.h"
#include <math.h>


Read_Sam_For_Training::Read_Sam_For_Training(const string & _sam_file_name,const string & _contig_file_name)
	:sam_file_name(_sam_file_name),contig_file_name(_contig_file_name)
{
	sam_info_1.resize(SIZE_OF_SAM);
	sam_info_2.resize(SIZE_OF_SAM);
}

Read_Sam_For_Training::Read_Sam_For_Training()
{
	throw runtime_error("Input sam file is empty!");
	exit(EXIT_FAILURE);
}

Read_Sam_For_Training::~Read_Sam_For_Training()
{
	//to do
}

string Read_Sam_For_Training::reverse(const string &s)
{
	string ret;
	if(s.empty())return ret;
	for(int i  = s.size()-1; i >= 0 ; --i)
	{
		switch(s[i])
		{
			case 'A':
				ret += 'T';
				break;
			case 'T':
				ret += 'A';
				break;
			case 'C':
				ret += 'G';
				break;
			case 'G':
				ret += 'C';
				break;
			default:
				ret += 'N';
				break;
		}
	}
	assert(ret.size() == s.size());
	return ret;
}


void Read_Sam_For_Training::update_insert(long insert_len)
{
	
	insert_len = abs(insert_len);
	
	//insert_len >= 0
	if (insert_len <= count_struct.max_insert_len)
	{
		count_struct.insert_lens[insert_len]++;
	}
	else
	{
		if (insert_len > MAX_INSERT_LEN)
		{
			count_struct.unused_read_count++;
			return;
		}

	}
	return;
}

void Read_Sam_For_Training::update_error_type(vector<string> &sam_info)
{

	string read = sam_info.at(READ_SEQ);
	int cur_read_len = read.size();

	for( int i = 0 ; i < cur_read_len ; ++i )
	{
		switch(read[i])
		{
			case 'A':
				count_struct.base_counts[A]++;break;

			case 'T':
				count_struct.base_counts[T]++;break;

			case 'C':
				count_struct.base_counts[C]++;break;

			case 'G':
				count_struct.base_counts[G]++;break;

			default:
				count_struct.base_counts[N]++;break;
		}
	}

	++(count_struct.read_lens[cur_read_len - 1]);

	string no_error_MD(int2Str(cur_read_len));

	string ref_qual = sam_info[MD_CHAR];

	if (ref_qual == no_error_MD)
	{
		return;
	}
	else
	{
		++(count_struct.error_read_count);
	}

	int ref_qual_len = ref_qual.size();
	//bool f_or_r = ((atoi(sam_info[LABEL].c_str()) >> 4) & 1)?true:false;
	bool f_or_r = false;

	string read_qual = sam_info[CIGAR];

	//ref_pos and read_pos is the cur pos of ref_seq and read
	int ref_pos = 0; // record the map pos 
	int read_pos = 0;

	int start = 0;

	char map_type;
	int cur_type_base_count = 0;

	//record inserts of each pos in read 
	vector<int> inserts_in_each_pos(cur_read_len, 0 );
	int tail = read_qual.find_first_of(DELIM_OF_READ);

	while(tail != string::npos)
	{
		cur_type_base_count = atoi(read_qual.substr(start,tail-start).c_str());
		map_type = read_qual[tail];

		if(map_type == 'M')
		{
			read_pos += cur_type_base_count;
			ref_pos += cur_type_base_count;
		}
		else if (map_type == 'I' || map_type == 'S')
		{
			if (f_or_r)
			{
				++(count_struct.ins_pos[ cur_read_len - read_pos - 1 ]);
				++(count_struct.ins_len[ cur_type_base_count - 1 ]);
			}
			else
			{
				++(count_struct.ins_pos[read_pos]);
				++(count_struct.ins_len[ cur_type_base_count ]);

			}
			inserts_in_each_pos[read_pos] = cur_type_base_count;
			read_pos += cur_type_base_count;
		}
		else if (map_type == 'D')
		{
			//reverse
			if(f_or_r)
			{
				++(count_struct.del_pos[cur_read_len - read_pos + 1]);
				++(count_struct.del_len[cur_type_base_count]);
			}
			else
			{
				++(count_struct.del_pos[read_pos]);
				++(count_struct.del_len[cur_type_base_count]);

			}
		}
		else
		{
			//nothing
		}
		start = tail+1;
		tail = read_qual.find_first_of(DELIM_OF_READ, start);

	}

	//quality
	cur_type_base_count = 0;
	ref_pos = 0;
	read_pos = 0;
	start = 0;
	tail = ref_qual.find_first_of(DELIM_OF_REF);
	

	while(tail != string::npos)
	{
		cur_type_base_count = atoi(ref_qual.substr(start,tail-start).c_str());
		map_type = ref_qual[tail];

		if (map_type == '^')
		{
			for(int i = tail; i < ref_qual_len ; ++i)
			{
				if (ref_qual[i] == 'A' ||
						ref_qual[i] == 'T' ||
						ref_qual[i] == 'C' ||
						ref_qual[i] == 'G' || 
						ref_qual[i] == 'N')
				{
					++tail;
				}
				else
				{
					break;
				}
			}
		}
		else if (map_type == 'A' || 
				map_type == 'T' ||
				map_type == 'C' ||
				map_type == 'G' ||
				map_type == 'N')
		{
			read_pos += cur_type_base_count;

			int cur_inserts_count = accumulate(inserts_in_each_pos.begin(), 
					inserts_in_each_pos.begin() + read_pos, 0);
			
			char to = read[cur_inserts_count + read_pos];

			//reverse
			if(f_or_r)
			{	
				++(count_struct.error_pos[cur_read_len - read_pos - 1]);
			}
			else
			{
				++(count_struct.error_pos[read_pos]);
			}
			int f, t;
			switch(map_type)
			{
				case 'A': f = A; break;
				case 'T': f = T; break;
				case 'C': f = C; break;
				case 'G': f = G; break;
				default:  f = N; break;
			}
			switch(to)
			{
				case 'A': t = A; break;
				case 'T': t = T; break;
				case 'C': t = C; break;
				case 'G': t = G; break;
				default:  t = N; break;
			}

			if(f == t)
			{
				//do nothing
			}
			else
			{
				++(count_struct.error_types[f][t]);
			}
			++read_pos;
		}

		start = tail+1;
		tail = ref_qual.find_first_of(DELIM_OF_REF, start);
	}

	return;
}

void Read_Sam_For_Training::read_contig_file()
{
	//ifstream contig_file_desc((TEMP_FILE_DIR + contig_file_name).c_str());
	ifstream contig_file_desc(contig_file_name.c_str());
	if (!contig_file_desc)
	{
		throw runtime_error("No such file or directory!!!");
		exit(1);
	}	
	string line;

	while(getline(contig_file_desc,line))
	{
		if (line[0] == '>')
		{
			continue;
		}
		else
		{
			(count_struct.contigs).push_back(line);
			(count_struct.contig_lens).push_back(line.size());
			(count_struct.total_contig_len) += line.size();
		}

	}
	contig_file_desc.close();
	contig_file_desc.clear();

	cout << "contigs.size() = " << count_struct.contigs.size() << endl;
	cout << "contig_lens.size() = " << count_struct.contig_lens.size() << endl;
	cout << count_struct.contig_lens[0] << "\t" << count_struct.contig_lens[count_struct.contig_lens.size()-1] << endl;
	cout << "total_contig_len = " << count_struct.total_contig_len << endl;

}

void Read_Sam_For_Training::find_max_common_seq1()
{

	string read = sam_info_1.at(READ_SEQ);
	string seq_name = sam_info_1.at(SEQ_NAME);

	vector<string> fields = (Split_String(seq_name)).split('_');
	int seq_index;

	seq_index = atoi(fields[1].c_str());

	long start_pos = atoi(sam_info_1.at(LEFT_END_POS).c_str());
	int read_len = sam_info_1.at(READ_SEQ).size();


	assert(count_struct.contigs.size() > seq_index);

cout << count_struct.contigs[seq_index].length() <<"-"<<start_pos-1 <<"-"<< read_len+DELTA << endl;
	string ref = count_struct.contigs[seq_index].substr(start_pos-1, read_len+ DELTA);


	//swap(read,ref);

	int row = read.size() + 1;
	int column = ref.size() + 1;

	//vector<vector<int> > matrix(row, vector<int>(column,0));
	//vector<vector<int> > back_track(row, vector<int>(column,-1));
	
	matrix.resize(row, vector<int>(column, -1));
	back_track.resize(row, vector<int>(column, -1));
	
	// initialization
	for(int j=0; j<column; ++j)
		matrix[0][j] = -3 * j;

	for(int i=0; i<row; ++i)
		matrix[i][0] = -3 * i;

	/*
	for(int i = 0; i < row-1 ; ++i)
	{
		for(int j = 0; j < column-1; ++j)
		{
			if (read[i] == ref[j])
			{
				matrix[i+1][j+1] = matrix[i][j] + 1;
				back_track[i+1][j+1] = 1;
			}
			if (matrix[i+1][j+1] < matrix[i+1][j])
			{
				matrix[i+1][j+1] = matrix[i+1][j];
				back_track[i+1][j+1] = 2;
			}
			if (matrix[i+1][j+1] < matrix[i][j+1])
			{
				matrix[i+1][j+1] = matrix[i][j+1];
				back_track[i+1][j+1] = 0;
			}
		}
	}
	*/
	// matrix[i][j] --- s[i] t[j] maximum matched count
	for(int i = 0; i < row-1 ; ++i)
	{
		for(int j = 0; j < column-1; ++j)
		{
			// match bonus = 1
			if (read[i] == ref[j])
			{
				matrix[i+1][j+1] = matrix[i][j]+1;
				back_track[i+1][j+1] = 1;
			}
			else { //mismatch bonus = 0
				matrix[i+1][j+1] = matrix[i][j]-2;
				back_track[i+1][j+1] = 1;
			}

			// bonus = -3 
			if (matrix[i+1][j+1] < matrix[i+1][j] - 3)
			{
				matrix[i+1][j+1] = matrix[i+1][j] - 3;
				back_track[i+1][j+1] = 2;
			}

			if (matrix[i+1][j+1] < matrix[i][j+1] - 3)
			{
				matrix[i+1][j+1] = matrix[i][j+1] - 3;
				back_track[i+1][j+1] = 0;
			}
		}
	}


	int i = row-1;
	int j = column-1;

	vector<int> path(read.size()+2,-1);
	//path.resize(column + 1,-1);

	while(back_track[i][j] != -1)
	{
		if (back_track[i][j] == 0)
		{
			--i;
		}
		else if (back_track[i][j] == 1)
		{
			path[i] = j;
			--i;--j;
		}
		else
		{
			--j;
		}
	}
	
	int left = 0, right = 0;
	int minus_one_count = 0, match_count = 0;
	string cigar;

	//boundary treatment
	path[0] = 0;
	path[path.size()-1] = path.size() - 1;

	i = 1;
	match_count = 0;
	while(i < path.size() - 1)
	{
		// increase one each time
		while(i < path.size() - 2 && path[i] != -1 && path[i+1] == path[i] + 1)
		{
			++ match_count;
			++ i;
		}
		if(i == path.size() - 2)
		{
			++ match_count;
			cigar += int2Str(match_count) + 'M';
			break;
		}
		else if(path[i] != -1 &&  path[i+1] != -1 && path[i+1] != path[i] + 1)
		{
			++ match_count;
			cigar += int2Str(match_count) + 'M';
			match_count = 0;

			cigar += int2Str(path[i+1] - path[i] - 1) + 'D';
			++ i;
			continue;
		}
		else if(path[i] != -1 && path[i+1] == -1)
		{
			++ match_count;
			++ i;
		}

		// [left+1, right-1] = {-1, -1 ,,,,-1}
		left = i-1;
		while(++i < path.size() -1  && path[i] == -1);
		right = i;

		int delta = (path[right] - path[left] - 1) - (right-left-1);
		// match
		if(delta == 0)
		{
			//cigar += int2Str(match_count + delta) + 'M';
			++ match_count;
		}
		// mismatch + delete
		else if(delta > 0)
		{

			if(right-left+1 + match_count > 0)
				cigar += int2Str(right-left-1 + match_count) + 'M';

			if(right != path.size() - 1)
				cigar += int2Str(delta) + 'D';

			match_count = 0;
		}
		// mismatch + insert
		else
		{

			if(path[right]-path[left]-1+match_count > 0)
				cigar += int2Str(path[right]-path[left]-1 + match_count) + 'M';

			if(right != path.size() - 1)
				cigar += int2Str(-delta) + 'I';

			match_count = 0;
		}

	}

	/*
	cout << "path : " << endl;
	copy(path.begin(),path.end(),ostream_iterator<int>(cout, "\t"));
	cout << endl;
	*/

	//from cigar and read seq to get the reference quality
	string quality;
	int tail = cigar.find_first_of(DELIM_OF_READ);

	//find quality
	int cur_type_base_count = 0;
	char map_type;
	int count = 0;
	int start = 0;

	// index of string s and t
	i = 0; j = 0;

	while(tail != string::npos)
	{
		map_type = cigar[tail]; 

cout << cigar.length() <<" "<< start <<" " <<tail-start <<endl;
		cur_type_base_count = atoi(cigar.substr(start,tail-start).c_str());

		if (map_type == 'M')
		{
			for(int k = 0 ;k < cur_type_base_count;++k)
			{
				// match
				if (read[i] == ref[j])
				{
					++ count;
				}
				else
				{ // mismatch
					if(count > 0){
						quality += int2Str(count);
						count = 0;
					}
					quality += ref[j];
				}
				++i;++j;
			}	
			if(count > 0){
				quality += int2Str(count);
				count = 0;
			}

		}
		else if (map_type == 'D')
		{			
			if(count > 0){
				quality += int2Str(count);
				count = 0;
			}

			quality += '^';

			for(int k = 0; k<cur_type_base_count; ++k)
			{
				quality += ref[j];
				++j;
			}
		}
		else if (map_type == 'I')
		{
			i += cur_type_base_count;
		}
		else
		{
			//do nothing
		}

		start = tail+1;
		tail = cigar.find_first_of(DELIM_OF_READ, start);
	}
	if(count > 0){
		quality += int2Str(count);
		count = 0;
	}

	sam_info_1.at(CIGAR) = cigar;
	sam_info_1.at(MD_CHAR) = quality;

#ifdef DEBUG
	cout << read << endl;
	cout << ref << endl;
	cout << cigar << endl;
	cout << quality << endl;
#endif

}


void Read_Sam_For_Training::find_max_common_seq2()
{

	string read = sam_info_2[READ_SEQ];
	string seq_name = sam_info_2[SEQ_NAME];

	vector<string> fields = (Split_String(seq_name)).split('_');
	int seq_index;

	//DEBUG
	if (fields.size() == 1)
	{
		cout << "error!!!" << endl;
		seq_index = 0;
	}
	//DEBUG
	else
	{
		seq_index = atoi(fields[1].c_str());
	}

	long start_pos = atoi(sam_info_2[LEFT_END_POS].c_str());
	int read_len = sam_info_2[READ_SEQ].size();

	assert(count_struct.contigs.size() > seq_index);
	string ref = count_struct.contigs[seq_index].substr(start_pos-1, read_len+DELTA);

	int row = read.size() + 1;
	int column = ref.size() + 1;

	//vector<vector<int> > matrix(row, vector<int>(column,0));
	//vector<vector<int> > back_track(row, vector<int>(column,-1));
	
	matrix.resize(row,vector<int>(column, -1));
	back_track.resize(row,vector<int>(column, -1));
	
	// initialization
	for(int j=0; j<column; ++j)
		matrix[0][j] = -3 * j;

	for(int i=0; i<row; ++i)
		matrix[i][0] = -3 * i;

	// matrix[i][j] --- s[i] t[j] maximum matched count
	for(int i = 0; i < row-1 ; ++i)
	{
		for(int j = 0; j < column-1; ++j)
		{
			// match bonus = 1
			if (read[i] == ref[j])
			{
				matrix[i+1][j+1] = matrix[i][j]+1;
				back_track[i+1][j+1] = 1;
			}
			else { //mismatch bonus = 0
				matrix[i+1][j+1] = matrix[i][j]-2;
				back_track[i+1][j+1] = 1;
			}

			// bonus = -3 
			if (matrix[i+1][j+1] < matrix[i+1][j] - 3)
			{
				matrix[i+1][j+1] = matrix[i+1][j] - 3;
				back_track[i+1][j+1] = 2;
			}

			if (matrix[i+1][j+1] < matrix[i][j+1] - 3)
			{
				matrix[i+1][j+1] = matrix[i][j+1] - 3;
				back_track[i+1][j+1] = 0;
			}
		}
	}

	/*
	for(int i = 0; i < row-1 ; ++i)
	{
		for(int j = 0; j < column-1; ++j)
		{
			if (read[i] == ref[j])
			{
				matrix[i+1][j+1] = matrix[i][j] + 1;
				back_track[i+1][j+1] = 1;
			}
			if (matrix[i+1][j+1] < matrix[i+1][j])
			{
				matrix[i+1][j+1] = matrix[i+1][j];
				back_track[i+1][j+1] = 2;
			}
			if (matrix[i+1][j+1] < matrix[i][j+1])
			{
				matrix[i+1][j+1] = matrix[i][j+1];
				back_track[i+1][j+1] = 0;
			}
		}
	}
	*/

	int i = row-1;
	int j = column-1;

	vector<int> path(read.size()+2,-1);
	//path.resize(column + 1,-1);

	while(back_track[i][j] != -1)
	{
		if (back_track[i][j] == 0)
		{
			--i;
		}
		else if (back_track[i][j] == 1)
		{
			path[i] = j;
			--i;--j;
		}
		else
		{
			--j;
		}
	}
	
	int left = 0, right = 0;
	int minus_one_count = 0, match_count = 0;
	string cigar;

	//boundary treatment
	path[0] = 0;
	path[path.size()-1] = path.size() - 1;

	i = 1;
	match_count = 0;
	while(i < path.size() - 1)
	{
		// increase one each time
		while(i < path.size() - 2 && path[i] != -1 && path[i+1] == path[i] + 1)
		{
			++ match_count;
			++ i;
		}
		if(i == path.size() - 2)
		{
			++ match_count;
			cigar += int2Str(match_count) + 'M';
			break;
		}
		else if(path[i] != -1 &&  path[i+1] != -1 && path[i+1] != path[i] + 1)
		{
			++ match_count;
			cigar += int2Str(match_count) + 'M';
			match_count = 0;

			cigar += int2Str(path[i+1] - path[i] - 1) + 'D';
			++ i;
			continue;
		}
		else if(path[i] != -1 && path[i+1] == -1)
		{
			++ match_count;
			++ i;
		}

		// [left+1, right-1] = {-1, -1 ,,,,-1}
		left = i-1;
		while(++i < path.size() -1  && path[i] == -1);
		right = i;

		int delta = (path[right] - path[left] - 1) - (right-left-1);
		// match
		if(delta == 0)
		{
			//cigar += int2Str(match_count + delta) + 'M';
			++ match_count;
		}
		// mismatch + delete
		else if(delta > 0)
		{

			if(right-left+1 + match_count > 0)
				cigar += int2Str(right-left-1 + match_count) + 'M';

			if(right != path.size() - 1)
				cigar += int2Str(delta) + 'D';

			match_count = 0;
		}
		// mismatch + insert
		else
		{

			if(path[right]-path[left]-1+match_count > 0)
				cigar += int2Str(path[right]-path[left]-1 + match_count) + 'M';

			if(right != path.size() - 1)
				cigar += int2Str(-delta) + 'I';

			match_count = 0;
		}

	}


	//from cigar and read seq to get the reference quality
	string quality;
	int tail = cigar.find_first_of(DELIM_OF_READ);

	//find quality
	int cur_type_base_count = 0;
	char map_type;
	int count = 0;
	int start = 0;
	int k = 0;

	// index of string s and t
	i = 0; j = 0;

	while(tail != string::npos)
	{
		map_type = cigar[tail]; 
		cur_type_base_count = atoi(cigar.substr(start,tail-start).c_str());

		if (map_type == 'M')
		{
			for(int k=0; k < cur_type_base_count; ++k)
			{
				if (read[i] == ref[j])
				{
					++ count;
				}
				else
				{
					if(count > 0){
						quality += int2Str(count);
						count = 0;
					}
					quality += ref[j];
				}
				++i;++j;
			}
			
			if(count > 0){
				quality += int2Str(count);
				count = 0;
			}
		}
		else if (map_type == 'D')
		{		
			if(count > 0){
				quality += int2Str(count);
				count = 0;
			}

			quality += '^';

			for(k = 0; k<cur_type_base_count; ++k)
			{
				quality += ref[j];
				++j;
			}
		}
		else if (map_type == 'I')
		{
			i += cur_type_base_count;
		}
		else
		{
			//do nothing
		}

		start = tail+1;
		tail = cigar.find_first_of(DELIM_OF_READ, start);
	}
	if(count > 0){
		quality += int2Str(count);
		count = 0;
	}

	sam_info_2[CIGAR] = cigar;
	sam_info_2[MD_CHAR] = quality;

	matrix.clear();
	back_track.clear();

#ifdef DEBUG
	cout << read << endl;
	cout << ref << endl;
	cout << cigar << endl;
	cout << quality << endl;
#endif
}


void Read_Sam_For_Training::set_two_quality()
{

	find_max_common_seq1();
	find_max_common_seq2();
}

void Read_Sam_For_Training::add_missing_info()
{
	//cout << "add missing info..." << endl;

	set_two_quality();

	long left_pos_1 = atoi(sam_info_1.at(LEFT_END_POS).c_str());
	long left_pos_2 = atoi(sam_info_2.at(LEFT_END_POS).c_str());


	sam_info_1.at(RIGHT_END_POS) = int2Str(left_pos_2);
	sam_info_2.at(RIGHT_END_POS) = int2Str(left_pos_1);


	int read_len = sam_info_2[READ_SEQ].size();

	string temp = int2Str(abs(left_pos_1 - left_pos_2) + read_len);
	sam_info_1.at(INSERT_LEN) = temp;
	sam_info_1.at(INSERT_LEN) = temp;

	sam_info_1.at(R_NEXT_INDEX) = R_NEXT; 
	sam_info_2.at(R_NEXT_INDEX) = R_NEXT; 
	
	long insert = atoi(sam_info_1.at(INSERT_LEN).c_str());
	if( insert > count_struct.max_insert_len)
	{
		count_struct.max_insert_len = insert;
	}
}


void Read_Sam_For_Training::filter_sam_file( const string &sam_file_name, const string &filter_sam_file_name )
{
	cout << "begin filter sam file for training..." << endl;

	cout << "sam_file_name : " << sam_file_name << endl;

	time_t start, end;
	time(&start);

	ifstream in(sam_file_name.c_str());
	if (!in)
	{
		throw runtime_error("No such file or directory!!");
		exit(1);
	}

	ofstream out(filter_sam_file_name.c_str());
	if (!out)
	{
		throw runtime_error("No such file or directory!!!");
		exit(1);
	}

	string line1,line2;

	int count = 0;
	size_t left, right;

	while(getline(in,line2))
	{
		++count;
		if (line2[0]  == '@' || line2[0] == ';') 
		{
			continue;
		}

		if (count % 1000000 == 0)
		{
			cout << count << endl;
		}

		if (line1.empty())
		{
			line1 = line2;

			right = line1.find_first_of("\t\n ");

			left = 0;
			int index = 0;
			while(string::npos != right)
			{
				sam_info_1.at(index) = line1.substr(left, right-left);
				++index;
				left = right + 1;
				right = line1.find_first_of("\t\n ",left);
			}
			sam_info_1.at(index) = line1.substr(left,right-left);
			int len1 = sam_info_1.at(READ_SEQ).size();
			//cout << len1 << endl;
#ifdef DEBUG
			assert(len1 == 36);
#endif
			int N_count = 0;
			for(int i = 0 ; i < len1 ; ++i)
			{
				if (sam_info_1.at(READ_SEQ).at(i) == 'N' || sam_info_1.at(READ_SEQ).at(i) == 'n')
				{
					N_count++;
				}
			}

			if (static_cast<double>(N_count/len1 < MAX_N_RATIO))
			{
				++(count_struct.total_read_count);
			}
			continue;
		}

		right = line2.find_first_of("\t\n ");
		left = 0;
		int	index = 0;
		while(string::npos != right)
		{
			sam_info_2[index++] = line2.substr(left, right-left);
			left = right + 1;
			right = line2.find_first_of("\t\n ",left);
		}
		sam_info_2[index] = line2.substr(left, right-left);
			
		//not a same pair
		if ( sam_info_1.at(READ_NAME) != sam_info_2.at(READ_NAME) || sam_info_1.at(SEQ_NAME) != sam_info_2.at(SEQ_NAME))
		{
			int len2 = sam_info_2.at(READ_SEQ).size();
			int N_count = 0;
			for(int i = 0 ; i < len2 ; ++i)
			{
				if (sam_info_2.at(READ_SEQ).at(i) == 'N')
				{
					N_count++;
				}
			}
			if(len2 > count_struct.max_read_len)
			{
				count_struct.max_read_len = len2;
			}

			if (static_cast<double>(N_count)/len2 < MAX_N_RATIO)
			{
				++(count_struct.unmapped_read_count);
				++(count_struct.total_read_count);
			}
			line1 = line2;
			sam_info_2[MD_CHAR] = "";
			sam_info_1.swap(sam_info_2);
			continue;
		}
		else
		{
			//a pair
			int t_size = sam_info_2.at(NM_CHAR).size();
			//multi map pair
			if ( sam_info_2.at(NM_CHAR).at(t_size - 1) > char('0' + 2))	
			{
				int len2 = sam_info_2.at(READ_SEQ).size();

				int N_count = 0;
				for(int i = 0 ; i < len2 ; ++i)
				{
					if (sam_info_2.at(READ_SEQ).at(i) == 'N')
					{
						N_count++;
					}
				}
				if(len2 > count_struct.max_read_len)
				{
					count_struct.max_read_len = len2;
				}

				if (static_cast<double>(N_count)/len2 < MAX_N_RATIO)
				{
					++(count_struct.unmapped_read_count);
					++(count_struct.total_read_count);
				}
				line1.clear();
				line2.clear();

				continue;
			}
			// a pair 

			if ( sam_info_1.at(R_NEXT_INDEX) == R_NEXT)
			{
				++(count_struct.total_read_count);

				line1.clear();
				line2.clear();
				sam_info_1.at(MD_CHAR) = int2Str(sam_info_1.at(READ_SEQ).size());
				sam_info_2[MD_CHAR] = int2Str(sam_info_2.at(READ_SEQ).size());

				int  i = 0;
				for( ; i < sam_info_1.size()-1 ; ++i)
				{
					out << sam_info_1.at(i) << "\t";
				}
				out << sam_info_1.at(i) << endl;

				i = 0;
				for( ; i < sam_info_2.size()-1 ; ++i)
				{
					out << sam_info_2.at(i) << "\t";
				}
				out << sam_info_2.at(i) << endl;

				continue;
			}
			//a pair but no more info
			add_missing_info();
			int  i = 0;
			for(; i < sam_info_1.size()-1;++i)
			{
				out << sam_info_1.at(i) << "\t";
			}
			out << sam_info_1.at(i) << endl;

			i = 0 ;
			for( ; i < sam_info_2.size()-1 ; ++i)
			{
				out << sam_info_2.at(i) << "\t";
			}
			out << sam_info_2.at(i) << endl;

			line1.clear();
			line2.clear();

			++(count_struct.total_read_count);
		}
	}

	in.close();
	in.clear();

	out.close();
	out.clear();

	time(&end);
	double seconds = difftime(end, start);
	cout << "filter file time = " <<  seconds << "  seconds " << endl;
}


void Read_Sam_For_Training::init_insert_counts(long _max_insert_len)
{
	//no need to verify
	
	if(_max_insert_len > MAX_INSERT_LEN)
	{
		count_struct.max_insert_len = MAX_INSERT_LEN;
		count_struct.insert_lens.resize(count_struct.max_insert_len + 1,1);
	}
	else
	{
		count_struct.max_insert_len = MAX_INSERT_LEN;
		count_struct.insert_lens.resize(count_struct.max_insert_len + 1,1);
	}
	
}

void Read_Sam_For_Training::init_error_types(int _max_read_len)
{
	count_struct.max_read_len = _max_read_len;

	count_struct.error_pos.resize(_max_read_len, 1);
	count_struct.ins_pos.resize(_max_read_len, 1);
	count_struct.ins_len.resize(_max_read_len, 1);
	count_struct.del_pos.resize(_max_read_len, 1);
	count_struct.del_len.resize(_max_read_len, 1);
	count_struct.read_lens.resize(_max_read_len, 0);
}



void Read_Sam_For_Training::update_parameter(Para &para)
{
	cout << "enter update..." << sam_file_name << endl;
	cout << "contig file name = " << contig_file_name << endl;

	read_contig_file();

	cout << "read contig file end..." << endl;

	//string filter_sam_file_name = TEMP_FILE_DIR + "filter_sam_for_training";
	string filter_sam_file_name = "filter_sam_for_training";

	cout << "begin filter sam file..." << endl;

	filter_sam_file(sam_file_name, filter_sam_file_name);

	cout << "end filter sam file..." << endl;

	//index form 1 not 0
	//++count_struct.max_insert_len;
	//	++count_struct.max_read_len;



	init_insert_counts(count_struct.max_insert_len);
	init_error_types(count_struct.max_read_len);


	ifstream sam_file_desc(filter_sam_file_name.c_str());
	if (!sam_file_desc)
	{
		throw runtime_error("No such file or directory!!!");
		exit(1);
	}


	bool f_or_r = false;

	string pre_read_name1("*");
	string pre_read_name2("*");

	string default_seq_name1("*");
	string default_seq_name2("*");

	string line1,line2;

	size_t left, right;

	string read_name1, read_name2;

	int count = 0; 
	vector<string> reads1, reads2;

	while(getline(sam_file_desc,line1))
	{	
		if(line1[0] == '@' || line1[0] == ';')
		{
			cout << "ERROR!!!!" << endl;
		}
		getline(sam_file_desc,line2);

		++count;
	
		//multi-mapped
		left = 0;
		right = line1.find_first_of("\t\n ");
		read_name1 = line1.substr(left, right-left);

		right = line2.find_first_of("\t\n ");
		read_name2 = line2.substr(left, right-left);
		
		assert(read_name1 == read_name2);

		if (read_name1 == pre_read_name1 && read_name2 == pre_read_name2)
		{

			reads1.push_back(line1);
			reads2.push_back(line2);

			continue;
		}
		else
		{
			assert(read_name1 == read_name2);
			//unique mapped
			if (reads1.size() == 1 && reads2.size() == 1)
			{
				
				string t_line1 = reads1[0];
				string t_line2 = reads2[0];

				left = 0;
				right = t_line1.find_first_of("\t\n ");
				int	index = 0;
				while(string::npos != right)
				{

					sam_info_1[index++] = t_line1.substr(left, right-left);
					left = right + 1;
					right = t_line1.find_first_of("\t\n ",left);
				}
				sam_info_1[index] = t_line1.substr(left,right-left);

				left = 0;
				right = t_line2.find_first_of("\t\n ");
				index = 0;
				while(string::npos != right)
				{
					sam_info_2[index++] = t_line2.substr(left,right-left);
					left = right+1;
					right = t_line2.find_first_of("\t\n ",left);
				}
				sam_info_2[index] = t_line2.substr(left,right-left);


				if (sam_info_1[MD_CHAR][0] != '^')
				{
					update_insert(atol(sam_info_1[INSERT_LEN].c_str()));
					count_struct.unique_mapped_count++;
				}
				if (sam_info_2[MD_CHAR][0] != '^')
				{
					update_insert(atol(sam_info_2[INSERT_LEN].c_str()));
					count_struct.unique_mapped_count++;
				}
				update_error_type(sam_info_1);
				update_error_type(sam_info_2);


				reads1.clear();
				reads2.clear();
			}

			reads1.push_back(line1);
			reads2.push_back(line2);
			pre_read_name1 =  read_name1;
			pre_read_name2 =  read_name2;
		}

	}
	sam_info_1.clear();
	sam_info_2.clear();

	sam_file_desc.close();
	sam_file_desc.clear();

	set_parameter(para);
}



void Read_Sam_For_Training::set_parameter(Para &para)
{
	//ofstream outfile((TEMP_FILE_DIR + "./count_struct").c_str());
	ofstream outfile("count_struct");

	outfile <<  count_struct;
	outfile << endl;

	outfile.close();
	outfile.clear();

	para.max_insert_len = count_struct.max_insert_len;
	para.max_read_len = count_struct.max_read_len;
	para.unique_mapped_count = count_struct.unique_mapped_count;
	para.unused_read_count = count_struct.unused_read_count;
	para.error_read_count = count_struct.error_read_count;
	para.contigs = count_struct.contigs;
	para.contig_lens = count_struct.contig_lens;
	para.total_contig_len = count_struct.total_contig_len;


	long temp_base_count = 0;
	for(int i = 0 ; i < ALPHABET ; ++i)
	{
		temp_base_count = 0;
		for(int j = 0 ; j < ALPHABET; ++j)
		{
			temp_base_count += count_struct.error_types.at(i).at(j);
		}

		for(int j = 0 ; j < ALPHABET; ++j)
		{
			para.base_error_type.at(i).at(j) = static_cast<double>(count_struct.error_types.at(i).at(j))/temp_base_count;
		}

		if (count_struct.base_counts.at(i) == 0)
		{
			para.base_error_rate.at(i) = ISNAN;
			continue;
		}
		para.base_error_rate.at(i) =  static_cast<double>(temp_base_count)/(double)(count_struct.base_counts.at(i));
	}

	//normalization
	double base_error_rate_sum = 0 ;
	for(int i = 0; i < ALPHABET-1 ; ++i)
	{
		base_error_rate_sum += para.base_error_rate.at(i);
	}

	for(int i = 0; i < ALPHABET-1 ; ++i)
	{
		para.base_error_rate.at(i) = (ALPHABET-1) * para.base_error_rate.at(i)/base_error_rate_sum;
	}
	//N is error!!
	para.base_error_rate.at(ALPHABET-1) = 1;

	assert(para.max_read_len == count_struct.max_read_len);
	assert(count_struct.read_lens.size() == count_struct.max_read_len);

	for(int i = para.max_read_len-1 ; i > 0 ; --i)
	{
		count_struct.read_lens.at(i-1) = count_struct.read_lens.at(i) + count_struct.read_lens.at(i-1);
	}

	(para.error_pos_dist).resize(para.max_read_len, 0);

	for(int i = 0; i < para.max_read_len ; ++i)
	{
		para.error_pos_dist.at(i) = static_cast<double>(count_struct.error_pos.at(i))/(count_struct.read_lens.at(i));
	}

	(para.ins_pos_dist).resize(para.max_read_len,0);

	for(int i = 0; i < para.max_read_len ; ++i)
	{
		para.ins_pos_dist.at(i) = static_cast<double>(count_struct.ins_pos.at(i))/(count_struct.read_lens.at(i));
	}

	(para.ins_len_dist).resize(para.max_read_len , 0);
	

	long temp_ins_count = 0 ;

	for(int i = 0; i < para.max_read_len ; ++i)
	{
		temp_ins_count += count_struct.ins_len.at(i);
	}

	for(int i = 0; i < para.max_read_len ; ++i)
	{
		para.ins_len_dist.at(i) = static_cast<double>(count_struct.ins_len.at(i))/temp_ins_count;
	}

	(para.del_pos_dist).resize(para.max_read_len, 0 );

	for(int i = 0; i < para.max_read_len ; ++i)
	{
		para.del_pos_dist.at(i) = static_cast<double>(count_struct.del_pos.at(i))/(count_struct.read_lens.at(i));
	}

	(para.del_len_dist).resize(para.max_read_len, 0 );

	long temp_del_count = 0 ;

	for(int i = 0; i < para.max_read_len ; ++i)
	{
		temp_del_count += count_struct.del_len.at(i);
	}
	for(int i = 0; i < para.max_read_len ; ++i)
	{
		para.del_len_dist.at(i) =  static_cast<double>(count_struct.del_len.at(i))/temp_del_count;
	}

	(para.insert_len_dist).resize(para.max_insert_len , 0);

	long int temp_insert_count = count_struct.unused_read_count;

	long temp_insert_sum = 0;
	for(int i = 0 ;i < para.max_insert_len ; ++i)
	{
		temp_insert_count += (count_struct.insert_lens.at(i));
		temp_insert_sum += (i*(count_struct.insert_lens.at(i)));
	}
	double insert_len_mean = static_cast<double>(temp_insert_sum)/temp_insert_count;

	temp_insert_sum = 0;
	for(int i = 0 ;i < para.max_insert_len ; ++i)
	{
		para.insert_len_dist.at(i) = static_cast<double>(count_struct.insert_lens.at(i))/temp_insert_count;

		temp_insert_sum += (count_struct.insert_lens.at(i))*(insert_len_mean-i) * (insert_len_mean-i);
	}

	double insert_len_var = temp_insert_sum/temp_insert_count; 


	(para.no_error_prob).resize(para.max_read_len , 0);

	double temp_no_error_prob = 1 ;
	for(int i = 0; i < para.max_read_len; ++i)
	{
		temp_no_error_prob *= (1 - para.error_pos_dist.at(i) - para.ins_pos_dist.at(i) - para.del_pos_dist.at(i));
		para.no_error_prob.at(i) = temp_no_error_prob;
	}

	cout << "max_insert_len = " << para.max_insert_len << endl;

	(para.effective_lens).resize( para.max_insert_len, -1);

	cout << "total_contig_len = " << count_struct.total_contig_len << endl;

	para.effective_lens.at(0) = count_struct.total_contig_len;

	assert(para.base_error_type.size() == ALPHABET);
	assert(para.base_error_rate.size() == ALPHABET);
	assert(para.error_pos_dist.size() == para.max_read_len);
	assert(para.ins_pos_dist.size() == para.max_read_len);
	assert(para.ins_len_dist.size() == para.max_read_len);
	assert(para.del_pos_dist.size() == para.max_read_len);
	assert(para.del_len_dist.size() == para.max_read_len);
	assert(para.insert_len_dist.size() == para.max_insert_len);
	assert(para.no_error_prob.size() == para.max_read_len);
	assert(para.effective_lens.size() == para.max_insert_len);
}

