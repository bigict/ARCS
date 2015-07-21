#include "Read_Sam_For_Test.h"
#include "Func.h"

ostream &operator<<(ostream &os, const Read_Sam_For_Test &obj)
{
	for(int i = 0 ; i < obj.gap_candidates_likelihood.size() ; ++i)
	{
		for( int j=0; j < obj.gap_candidates_likelihood[i].size() ; ++j)
		{
			os<< i << " " << j << " "  << obj.gap_candidates_likelihood[i][j] << " " << obj.gap_candidates_likelihood[i][j]/obj.read_mapped_count[i][j] << endl;
		}
		os << endl;
	}

	return os;
}


Read_Sam_For_Test::Read_Sam_For_Test()
{
	throw runtime_error("Input sam file is empty!");
	exit(EXIT_FAILURE);
}


Read_Sam_For_Test::Read_Sam_For_Test(const string& _gap_count_file_name, const string& _sam_file_name, const string &_contig_file_name)
:sam_file_name(_sam_file_name), contig_file_name(_contig_file_name)
{
	

	ifstream infile(_gap_count_file_name.c_str());
	if (!infile)
	{
		cout << _gap_count_file_name << endl;
		throw runtime_error("No such file or directory");
		exit(EXIT_FAILURE);
	}
/*
	if (!infile.eof())
	{
		infile >> gap_count;
	}
*/

	while(!infile.eof())
	{
		int temp;
		infile >> temp;
		if(infile.eof())
		{
			break;
		}
		gap_candi_count.push_back(temp);
	}
	
	infile.close();
	infile.clear();

	gap_count = gap_candi_count.size();
	assert(gap_count == gap_candi_count.size());

	cout << gap_count << "\t" << gap_candi_count.size() << endl;

	for(int i = 0; i < gap_count ; ++i)
	{

		vector<double> candidates(gap_candi_count[i],0);
		vector<long> read_mapped(gap_candi_count[i],0);
		
		gap_candidates_likelihood.push_back(candidates);
		read_mapped_count.push_back(read_mapped);

		all_contigs.push_back(vector<string>(gap_candi_count[i], ""));
	}

	assert(gap_candidates_likelihood.size() == gap_count);
	assert(read_mapped_count.size() == gap_count);

	sam_info_1.resize(SIZE_OF_SAM);
	sam_info_2.resize(SIZE_OF_SAM);
}




Read_Sam_For_Test::~Read_Sam_For_Test()
{
	//
}

//read scoring contig file
void Read_Sam_For_Test::read_contig_file()
{
	ifstream contig_file_desc(contig_file_name.c_str());
	if (!contig_file_desc)
	{
		throw runtime_error("No such file or directory!!!");
		exit(EXIT_FAILURE);
	}

	string line;
	string seq;

	int gap_index = -1, candi_index = -1;
	while(getline(contig_file_desc,line))
	{
		if (line[0] == ';' || line[0] == '@')
		{
			continue;
		}
		if (line[0] == '>')
		{
			if (!seq.empty())
			{
				assert(gap_index != -1 && candi_index != -1);

				all_contigs[gap_index][candi_index] = seq;
				seq.clear();
				gap_index = -1;
				candi_index = -1;
			}
			vector<string> fields = Split_String(line).split('_');
			gap_index = atoi(fields[1].c_str());
			candi_index = atoi(fields[2].c_str());
			continue;
		}
		else
		{
			seq += line;
		}
	}

	assert(!seq.empty());
	assert(gap_index != -1 && candi_index != -1);

	if (!seq.empty() && gap_index != -1 && candi_index != -1)
	{
		all_contigs[gap_index][candi_index] = seq;
	}

	contig_file_desc.close();
	contig_file_desc.clear();

	cout << "scoring contig number\t:\t " << all_contigs.size() << endl;
}

//update likelihood of currrent gap candidate
double Read_Sam_For_Test::compute_likelihood( Para &para,const int gap_index,const int candi_index)
{
	long insert_len = atoi(sam_info_1.at(INSERT_LEN).c_str());
	
	double insert_len_prob(0.0);
	if(insert_len == 0  || insert_len > MAX_INSERT_LEN || insert_len > para.max_insert_len)
	{
		insert_len_prob = 1/(para.unique_mapped_count);
	}
	else
	{
		insert_len_prob = para.insert_len_dist.at(insert_len);
	}

	double log_insert_len_prob(0.0);

	if( insert_len_prob < ISNAN || isnan( insert_len_prob ) )
	{
		insert_len_prob = ISNAN;
	} 

	log_insert_len_prob = log(insert_len_prob);
	
	double log_locate_prob = compute_locate_prob(insert_len , para, gap_index, candi_index);

	double log_err_prob1(0.0), log_err_prob2(0.0);

	log_err_prob1 = compute_error_prob1( para );

	log_err_prob2 = compute_error_prob2( para );

	double log_prob = log_locate_prob + log_insert_len_prob + log_err_prob1 + log_err_prob2;
	return log_prob;
}




double Read_Sam_For_Test::compute_error_prob1(const Para &para )
{
	
	string read = sam_info_1.at(READ_SEQ);
	string read_qual = sam_info_1.at(CIGAR);
	string ref_qual = sam_info_1.at(MD_CHAR);

	int cur_read_len = read.size();
	int read_qual_len = read_qual.size();
	int ref_qual_len = ref_qual.size();
    
	//bool f_or_r = ((atoi(sam_info_1[LABEL].c_str()) >> 4) & 1)?true:false;
	bool f_or_r = false;

	double log_error_prob = log(para.no_error_prob.at(cur_read_len-1));

	if (ref_qual[0] == '^')
	{
		return log_error_prob;
	}

	int ref_pos(0);
	int read_pos(0);
	int start(0);
	char map_type;
	int cur_type_base_count(0);

	vector<int> inserts_in_each_pos(cur_read_len,0);
	int tail = read_qual.find_first_of(DELIM_OF_READ);

	int index(0);

	while(tail != string::npos)
	{
		cur_type_base_count = atoi(read_qual.substr(start,tail-start).c_str());
		map_type = read_qual.at(tail);

		if(map_type == 'M')
		{
			read_pos += cur_type_base_count;
			ref_pos += cur_type_base_count;
		}
		else if (map_type == 'I' || map_type == 'S')
		{
			if(f_or_r)
			{
				index = cur_read_len - read_pos - 1;
			}
			else
			{
				index = read_pos;

			}
			log_error_prob = log_error_prob + log(para.ins_pos_dist.at(index)) + log(para.ins_len_dist.at(cur_type_base_count - 1)) - log(1 - (para.error_pos_dist.at(index)) - para.ins_pos_dist.at(index) - para.del_pos_dist.at(index));

			//current position has insertion
			inserts_in_each_pos.at(read_pos) = cur_type_base_count;
			read_pos += cur_type_base_count;
		}
		else if (map_type == 'D')
		{
			if(f_or_r)
			{
				index = cur_read_len - read_pos - 1;
			}
			else
			{
				index = read_pos;
			}	
			assert(para.del_pos_dist.size() == cur_read_len);
			assert(para.del_len_dist.size() == cur_read_len);
			assert(para.error_pos_dist.size() == cur_read_len);

			log_error_prob = log_error_prob + log(para.del_pos_dist.at(index)) + log(para.del_len_dist.at(cur_type_base_count - 1)) - log(1 - (para.error_pos_dist.at(index)) - para.ins_pos_dist.at(index) - para.del_pos_dist.at(index));

		}
		else
		{
			//nothing
		}
		start = tail + 1;
		tail = read_qual.find_first_of(DELIM_OF_READ, start);
	}

	cur_type_base_count = 0;
	ref_pos = 0;
	read_pos = 0;
	start = 0;
	tail = ref_qual.find_first_of(DELIM_OF_REF);
	f_or_r = false;

	while(tail != string::npos)
	{
		cur_type_base_count = atoi(ref_qual.substr(start, tail-start).c_str());
		map_type = ref_qual.at(tail);

		if(map_type == '^')
		{
			read_pos += cur_type_base_count;
			++tail;

			//from cur to find all ^ bases
			for(int  k = tail; k < ref_qual_len; ++k)
			{
				if (ref_qual[k] == 'A' ||
						ref_qual[k] == 'T' ||
						ref_qual[k] == 'C' ||
						ref_qual[k] == 'G' || 
						ref_qual[k] == 'N')
				{
					++tail;
				}
				else
				{
					--tail;
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

			int cur_inserts_count = accumulate(inserts_in_each_pos.begin(),	inserts_in_each_pos.begin() + read_pos , 0);

			char to = read.at(cur_inserts_count + read_pos);

			if(f_or_r == false)
			{
				index = read_pos;

			}
			else
			{
//				index = cur_read_len - read_pos - ref_pos;
			}

			log_error_prob = log_error_prob  + log(para.error_pos_dist.at(index)) - log( 1 - para.error_pos_dist.at(index) - para.ins_pos_dist.at(index) - para.del_pos_dist.at(index));

			int f,t;
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
				log_error_prob += (log( para.base_error_rate[f]) + log (para.base_error_type[f][t]));
			}
			read_pos += 1;
		}
		else
		{
			//break;	
			//do nothing
		}
		start = tail + 1;
		tail = ref_qual.find_first_of(DELIM_OF_REF, start);
	}

	return log_error_prob;
}


double Read_Sam_For_Test::compute_error_prob2(const Para &para )
{
	string read = sam_info_2.at(READ_SEQ);
	string read_qual = sam_info_2.at(CIGAR);
	string ref_qual = sam_info_2.at(MD_CHAR);

	int cur_read_len = read.size();
	int read_qual_len = read_qual.size();
	int ref_qual_len = ref_qual.size();

	//bool f_or_r = ((atoi(sam_info_2[LABEL].c_str()) >> 4) & 1)?true:false;
	bool f_or_r = false;

	double log_error_prob = log(para.no_error_prob.at(cur_read_len-1));

	if (ref_qual[0] == '^')
	{
		return log_error_prob;
	}

	int ref_pos(0);
	int read_pos(0);
	int start(0);
	char map_type;
	int cur_type_base_count(0);

	vector<int> inserts_in_each_pos(cur_read_len,0);
	int tail = read_qual.find_first_of(DELIM_OF_READ);

	int index(0);

	while(tail != string::npos)
	{
		cur_type_base_count = atoi(read_qual.substr(start,tail-start).c_str());
		map_type = read_qual.at(tail);

		if(map_type == 'M')
		{
			read_pos += cur_type_base_count;
			ref_pos += cur_type_base_count;
		}
		else if (map_type == 'I' || map_type == 'S')
		{
			//reverse
			if(f_or_r)
			{
				index = cur_read_len - read_pos - 1;
			}
			else
			{
				index = read_pos;
			}

			log_error_prob = log_error_prob + log(para.ins_pos_dist.at(index))	+ log(para.ins_len_dist.at(cur_type_base_count - 1))- log(1 - (para.error_pos_dist.at(index)) - para.ins_pos_dist.at(index) - para.del_pos_dist.at(index));
			
			inserts_in_each_pos.at(read_pos) = cur_type_base_count;

			read_pos += cur_type_base_count;
		}
		else if (map_type == 'D')
		{
			//reverse
			if(f_or_r)
			{
				index = cur_read_len - read_pos - 1;
			}
			else
			{
				index = read_pos ;
			}	
			log_error_prob = log_error_prob + log(para.del_pos_dist.at(index)) + log(para.del_len_dist.at(cur_type_base_count - 1)) - log(1 - (para.error_pos_dist.at(index)) - para.ins_pos_dist.at(index) - para.del_pos_dist.at(index));

		}
		else
		{
			//nothing
		}
		start = tail + 1;
		tail = read_qual.find_first_of(DELIM_OF_READ, start);
	}

	cur_type_base_count = 0;
	ref_pos = 0;
	read_pos = 0;
	start = 0;
	tail = ref_qual.find_first_of(DELIM_OF_REF);


	while(tail != string::npos)
	{
		cur_type_base_count = atoi(ref_qual.substr(start, tail-start).c_str());

		map_type = ref_qual.at(tail);

		if(map_type == '^')
		{
			//from cur to find all ^ bases
			read_pos += cur_type_base_count;
			++tail;

			for(int  k = tail; k < ref_qual_len; ++k)
			{
				if (ref_qual[k] == 'A' ||
						ref_qual[k] == 'T' ||
						ref_qual[k] == 'C' ||
						ref_qual[k] == 'G' || 
						ref_qual[k] == 'N')
				{
					++tail;
				}
				else
				{
					--tail;
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
			read_pos += cur_type_base_count ;
			//read_pos += cur_type_base_count;

			int cur_inserts_count = accumulate(inserts_in_each_pos.begin(),	inserts_in_each_pos.begin() + read_pos , 0);

			char to = read.at(cur_inserts_count + read_pos);
			
			int index = 0;
			if(f_or_r == false)
			{
				//index = read_pos - 1 - ref_pos;
				index = read_pos;
			}
			else
			{
				cout << "NEVER!!!" << endl;
				index = cur_read_len - read_pos - ref_pos;
			}

			log_error_prob = log_error_prob + log(para.error_pos_dist.at(index))- log(1-para.error_pos_dist.at(index)-para.ins_pos_dist.at(index)-para.del_pos_dist.at(index));
			int f,t;
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
				log_error_prob += (log( para.base_error_rate[f]) + log (para.base_error_type[f][t]));
			}
			read_pos += 1;
		}
		else
		{
			//break;	
			//do nothing
		}	
		start = tail + 1;
		tail = ref_qual.find_first_of(DELIM_OF_REF, start);
	}

	return log_error_prob;

}

//update the location site of proobability
double Read_Sam_For_Test::compute_locate_prob( const long insert_len, Para &para, const int gap_index, const int candi_index)
{
	double locate_prob(0);

	if(insert_len < 0)
	{
		locate_prob = ( 1/static_cast<double>(para.effective_lens[0]) );
	}

	long locate_site(0);
	if( insert_len > para.max_insert_len )
	{
		for(int i = 0 ; i < para.contig_lens.size() ; ++i)
		{
			locate_site += (para.contig_lens.at(i) - insert_len + 1 );
		}
		locate_prob = (1/static_cast<double>(locate_site));
	}
	else
	{
		locate_site = 0;

		for(int i = 0 ; i < para.contig_lens.size(); ++i)
		{
			locate_site += (para.contig_lens.at(i) - insert_len + 1);
		}
		locate_prob = (1/static_cast<double>(locate_site));
		(para.effective_lens).at(insert_len) = locate_site;
	}

	double log_locate_prob(0);

	if( locate_prob < ISNAN || isnan( locate_prob ) )
	{
		locate_prob = 1e-320;
	}
	log_locate_prob = log(locate_prob);

	return log_locate_prob;
}	


void Read_Sam_For_Test::find_max_common_seq1(const int gap_index, const int candi_index)
{
	string read = sam_info_1.at(READ_SEQ);

	long start_pos = atoi(sam_info_1.at(LEFT_END_POS).c_str());
	int read_len = sam_info_1.at(READ_SEQ).size();

	string ref = all_contigs.at(gap_index).at(candi_index).substr(start_pos-1, read_len+ DELTA);


	int row = read.size() + 1;
	int column = ref.size() + 1;

	matrix.resize(row,vector<int>(column, -1));
	back_track.resize(row,vector<int>(column, -1));

	//initialize
	for(int j=0; j<column; ++j)
		matrix[0][j] = -3 * j;

	for(int i=0; i<row; ++i)
		matrix[i][0] = -3 * i;
/*
	for(int i = 0; i < row-1 ; ++i)
	{
		for(int j = 0; j < column-1; ++j)
		{
			if (ref_seq[i] == read[j])
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
	// match +1 mismatch -2 del/ins -3
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
			else { //mismatch bonus = -2
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

	int left(0), right(0);
	int minus_one_count(0), match_count(0);
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

	// index of string read(s) and ref(t)
	i = 0; j = 0;

	while(tail != string::npos)
	{
		map_type = cigar.at(tail);
		cur_type_base_count = atoi(cigar.substr(start,tail-start).c_str());

		if (map_type == 'M')
		{
			for(int k = 0 ; k < cur_type_base_count ; ++k)
			{
				//match
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

			for(int k = 0; k < cur_type_base_count; ++k)
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
	cout << read << endl;
	cout << ref << endl;
	cout << cigar << endl;
	cout << quality << endl;

}


void Read_Sam_For_Test::find_max_common_seq2(int gap_index, int candi_index)
{

	string read = sam_info_2.at(READ_SEQ);

	long start_pos = atoi(sam_info_2.at(LEFT_END_POS).c_str());
	int read_len = sam_info_2.at(READ_SEQ).size();

	string ref = all_contigs.at(gap_index).at(candi_index).substr(start_pos-1, read_len+DELTA);


	int row = read.size() + 1;
	int column = ref.size() + 1;

	matrix.resize(row,vector<int>(column, -1));
	back_track.resize(row,vector<int>(column, -1));


	//initialize
	for(int j=0; j<column; ++j)
		matrix[0][j] = -3 * j;

	for(int i=0; i<row; ++i)
		matrix[i][0] = -3 * i;
/*

	for(int i = 0; i < row-1 ; ++i)
	{
		for(int j = 0; j < column-1; ++j)
		{
			if (ref_seq[i] == read[j])
			{
				matrix[i+1][j+1] = matrix[i][j] + 1;
				back_track[i+1][j+1] = 1;
			}
			if (matrix[i+1][j+1] < matrix[i+1][j])
			{
				matrix[i+1][j+1] = matrix[i+1][j];
				back_track[i+1][j+1] = 2; }
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
			else { //mismatch bonus = -2
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
			for( k = 0; k < cur_type_base_count ; ++k)
			{
				if (read[i] == ref[j])
				{
					++ count;
				}
				else
				{
					if(count > 0){
						quality += int2Str(count);
						//cout << count << endl;
						count = 0;
					}
					quality += ref[i];
				}
				++i;++j;
			}
			
		   if(count > 0){
			   quality += int2Str(count);
		//	   cout << count << endl;
			   count = 0;
		   }
		}
		else if (map_type == 'D')
		{		
			if(count > 0){
				quality += int2Str(count);
				//cout << count << endl;
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
		//cout << count << endl;
		count = 0;
	}

	sam_info_2.at(CIGAR) = cigar;
	sam_info_2.at(MD_CHAR) = quality;
	matrix.clear();
	back_track.clear();

#ifdef DEBUG
	cout << read << endl;
	cout << ref << endl;
	cout << cigar << endl;
	cout << quality << endl;
#endif
}


void Read_Sam_For_Test::set_two_quality(const int gap_index,const int candi_index)
{

	find_max_common_seq1(gap_index, candi_index);
	find_max_common_seq2(gap_index, candi_index);
}

void Read_Sam_For_Test::add_missing_info(const int gap_index, const int candi_index)
{

	set_two_quality(gap_index, candi_index);

	long left_pos_1 = atoi(sam_info_1.at(LEFT_END_POS).c_str());
	long left_pos_2 = atoi(sam_info_2.at(LEFT_END_POS).c_str());


	sam_info_1.at(RIGHT_END_POS) = int2Str(left_pos_2);
	sam_info_2.at(RIGHT_END_POS) = int2Str(left_pos_1);

	int read_len = sam_info_2.at(READ_SEQ).size();

	string temp = int2Str((abs( left_pos_1 - left_pos_2 ) + read_len));

	sam_info_1.at(INSERT_LEN) = temp;
	sam_info_2.at(INSERT_LEN) = temp;

	sam_info_1.at(R_NEXT_INDEX) = R_NEXT; 
	sam_info_2.at(R_NEXT_INDEX) = R_NEXT; 
}

//filter the sam file of scoring contig 
void Read_Sam_For_Test::filter_sam_file( const string& sam_file_name, const string& filter_sam_file_name )
{	
    cout << "Read the sam file of scoring contigs..." << endl;
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

	int count(0);
	size_t left, right;

	while(getline(in,line2))
	{
		++count;
		if (line2[0]  == '@' || line2[0] == ';') 
        {
          continue;
        }

		if (line1.empty())
		{
			line1 = line2;

			right = line1.find_first_of("\t");
			left = 0;
			int index = 0;
			while(string::npos != right)
			{
				sam_info_1.at(index) = line1.substr(left, right-left);
				++index;
				left = right + 1;
				right = line1.find_first_of("\t",left);
			}
			sam_info_1.at(index) = line1.substr(left,right-left);
			continue;
		}
		right = line2.find_first_of("\t");
		left = 0;
		int	index = 0;
		while(string::npos != right)
		{
			sam_info_2.at(index) = line2.substr(left, right-left);
			++index;
			left = right + 1;
			right = line2.find_first_of("\t", left);
		}
		sam_info_2.at(index)= line2.substr(left, right-left);

        //mapped paired reads
		if (sam_info_1.at(READ_NAME) == sam_info_2.at(READ_NAME) && 
            sam_info_1.at(SEQ_NAME) == sam_info_2.at(SEQ_NAME) && 
            sam_info_1.at(R_NEXT_INDEX) == R_NEXT && 
            sam_info_2.at(R_NEXT_INDEX ) == R_NEXT)
		{	
			vector<string> fields = (Split_String(sam_info_1.at(SEQ_NAME))).split('_');
			assert(fields.size() == 3);
			int gap_index = atoi(fields[1].c_str());
			int candi_index = atoi(fields[2].c_str());

            //set MD_char and cigar
			set_two_quality(gap_index, candi_index);
            
            //output filter result
			int i = 0 ;
			for(; i < sam_info_1.size()-1;++i)
			{
				out << sam_info_1.at(i) << "\t";
			}
			out << sam_info_1.at(i) << endl;
			i = 0 ;
			for(; i < sam_info_2.size()-1;++i)
			{
				out << sam_info_2.at(i) << "\t";
			}
			out << sam_info_2.at(i) << endl;

		}
		sam_info_1.swap(sam_info_2);
	}
	in.close();
	in.clear();

	out.close();
	out.clear();

	time(&end);
	double seconds = difftime(end, start);
	cout << "Filter sam file for scoring :  " <<  seconds << "  seconds " << endl;
}



//update corresponding Candi_Info struct according to gap_index and candi_index
void Read_Sam_For_Test::scoring( Para &para)
{
    cout << "Enter scoring..." << endl;
	read_contig_file();
	ifstream sam_file_desc(sam_file_name.c_str());
	if (!sam_file_desc)
	{
		cout << sam_file_name << "No such file or directory!!" << endl;
		exit(EXIT_FAILURE);
	}

    char buf[10];
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    string fileName;
    fileName += string(buf);
    fileName += string("mer");
    fileName += ".filter_sam_for_scoring";


	//string filter_sam_file_name = TEMP_FILE_DIR + "filter_sam_for_test";
	string filter_sam_file_name = fileName;
	filter_sam_file(sam_file_name , filter_sam_file_name);

	string line1,line2;
	size_t left, right;

	long gap_index(-1);
	long candi_index(-1);

	double log_prob(0);
	int count(0);	

	sam_file_desc.close();
	sam_file_desc.clear();

	sam_file_desc.open(filter_sam_file_name.c_str());

	if (!sam_file_desc)
	{
		cout << "open filter sam file error!!!" << endl;
		exit(1);
	}

    //read sam file and extract the useful file name
	while(getline(sam_file_desc, line1))
	{
		if(line1[0] == '@' || line1[0] == ';')
		{
			continue;
		}
		getline(sam_file_desc, line2);
        
        //extract information of sam
		left = 0;
		right = line1.find_first_of("\t");
		int	index = 0;
		while(string::npos != right)
		{
			sam_info_1[index] = line1.substr(left, right-left);
			++index;
			left = right + 1;
			right = line1.find_first_of("\t",left);
		}

		left = 0;
		right = line2.find_first_of("\t");
		index = 0;
		while(string::npos != right)
		{
			sam_info_2[index] = line2.substr(left, right-left);
			++index;
			left = right + 1;
			right = line2.find_first_of("\t",left);
		}

		vector<string> fields = (Split_String(sam_info_1[SEQ_NAME])).split('_');

		assert(fields.size() == 3);

		gap_index = atoi(fields[1].c_str());
		candi_index = atoi(fields[2].c_str());
        
        //compute likelihood
		double log_prob = compute_likelihood( para , gap_index, candi_index );
	
		gap_candidates_likelihood[gap_index][candi_index] += log_prob;
		read_mapped_count[gap_index][candi_index]++;
	}

	sam_file_desc.close();
	sam_file_desc.clear();
	cout << "begin output right choice..." << endl;
	output_final_scaffolds();


}

bool myFunc(double lhs, double rhs)
{
	return (lhs > rhs);
}



void Read_Sam_For_Test::output_final_scaffolds()
{
	//string right_gap_index_name = TEMP_FILE_DIR + "each_gap_right_index";
	string right_gap_index_name = "each_gap_right_index";
	ofstream out(right_gap_index_name.c_str());

	double max_index = 0 ;
	vector<double>::iterator iter;
	vector<int> right_choices;

	for(int i = 0 ; i < gap_count; ++i)
	{
		iter = max_element(
				gap_candidates_likelihood[i].begin(),
				gap_candidates_likelihood[i].end());
		max_index = iter - gap_candidates_likelihood[i].begin();
		out << i << " " << max_index << " " << *iter << endl;	
		right_choices.push_back(max_index);

	}
	assert(right_choices.size() == gap_count);

	out.close();
	out.clear();

	//string new_all_scaffolds_segments_file_name = TEMP_FILE_DIR + "new_all_scaffolds_segments_to_score";
	string new_all_scaffolds_segments_file_name = "new_all_scaffolds_segments_to_score";
	ifstream in(new_all_scaffolds_segments_file_name.c_str());
	string line;
	vector<vector<string> > scaffolds_segments;
	vector<string> seqs;

	while(getline(in, line))
	{
		if(line[0] == '>')
		{
			if(!seqs.empty())
			{
				scaffolds_segments.push_back(seqs);
			}
			seqs.clear();
		}
		seqs.push_back(line);
	}
	if (!seqs.empty())
	{
		scaffolds_segments.push_back(seqs);
	}
	cout << "scaffolds_segments.size() = " << scaffolds_segments.size()  << endl;

	in.close();
	in.clear();

	string seq;
	//string final_scaffold_file_name = TEMP_FILE_DIR + "final_scaffolds_seq";
	string final_scaffold_file_name = "final_scaffolds_seq";
	out.open(final_scaffold_file_name.c_str());

	int index = 0;
	int threshold = MU + var*3;

	for(int i  = 0; i < scaffolds_segments.size() ; ++i)
	{
		if (scaffolds_segments[i].size() == 1)
		{
			out << ">scaff " << scaffolds_segments[i].size() << endl; 
			out << scaffolds_segments[i][0] << endl;
			continue;
		}
		seq.clear();
		seq += scaffolds_segments[i][0];
		int  j = 1;
		for ( ; j  < scaffolds_segments[i].size() - 1 ; ++j)
		{
			seq += scaffolds_segments[i][j];
		//	seq.resize(seq.size() - threshold);
			seq += all_contigs[index][right_choices[index]];
		//	seq.resize(seq.size() - threshold);
			++index;
		}
		assert(j == scaffolds_segments[i].size()-1);
		seq += scaffolds_segments[i][j];
		out << ">scaff " << scaffolds_segments[i].size() << endl;
		out << seq << endl;
	}
	out.close();
	out.clear();
	assert(index == gap_count);

}

