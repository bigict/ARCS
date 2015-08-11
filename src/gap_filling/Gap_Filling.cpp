#include "Gap_Filling.h"
#include "Split_String.h"


#define TRAINING
#define RESULT

Gap_Filling::Gap_Filling(const string & _contig_file_name, const string &_scaffold_file_name, const string  &_initial_contig_file_name):edge_file_name(_contig_file_name),scaffold_file_name(_scaffold_file_name),initial_edge_file_name(_initial_contig_file_name)
{
	ifstream in(edge_file_name.c_str());
	if (!in) 
	{
		cerr << "[Info] " << edge_file_name << " No such file or directory!!!" << endl; 
		exit(1);
	}
	in.close();
	in.clear();

	in.open(scaffold_file_name.c_str());
	if (!in) 
	{
		cerr << "[Info] " << scaffold_file_name << " No such file or directory!!!" << endl; 
		exit(1);
	}
	in.close();
	in.clear();

	in.open(initial_edge_file_name.c_str());
	if (!in) 
	{
		cerr << "[Info] " << initial_edge_file_name << " No such file or directory!!!" << endl; 
		exit(1);
	}
	in.close();
	in.clear();

}

//dtor
Gap_Filling::~Gap_Filling()
{

}

//gap filling
void Gap_Filling::gap_filling()
{
	get_gap_info();
	uniq_candidate_gap_filling();
	multi_candidates_gap_filling();
	output_initial_scaffolds_seq();
	get_multi_candidates_gap_local_seq();
}


//Align two nerghboring contigs
int Gap_Filling::alignment(const string& suffix, const string & prefix)
{
	assert(prefix.size() == suffix.size());

	vector<int> score(prefix.size()+1,0);

	for(int i = 1 ; i < suffix.size() + 1 ;++i)
	{
		for(int j = prefix.size(); j > 0 ; --j)
		{
			if (prefix[j-1] == suffix[i-1])
			{
				score[j] = score[j-1]+1;
			}
			else
			{
				score[j] = 0;
			}
		}
	}

	int ret = 0;
	for(int i = 0; i < prefix.size()+1; ++i)
	{
		if (score[i] == i)
		{
			if (ret < score[i])
			{
				ret = score[i];
			}
		}
	}
	return ret;
}

static bool myFunc(const string &s1, const string &s2)
{
	return (s1.size() > s2.size());
}

class mycomparison				/*functor*/
{
    bool reverse;
public:
    mycomparison(const bool _reverse = false):reverse(_reverse){}

    bool operator() (const string& lhs, const string&rhs) const
    {
        if (reverse)
            return lhs.size() > rhs.size();
        else 
            return lhs.size() < rhs.size();
    }
};


//read contigs and contigs in condensed de Bruijn graph
void Gap_Filling::set_edges()
{
	ifstream infile(edge_file_name.c_str());
	if (!infile)
	{
		cerr << "[Info] " << edge_file_name << "No such file or directory!!" << endl;
		exit(EXIT_FAILURE);
	}

	string line;
	string seq;
	int cur_index(0);

	//get unique edges for training parameter
#ifdef TRAINING
    priority_queue<string, vector<string>, mycomparison> unique_contigs_for_training;
#endif

	vector<string> fields;
	int copy_number = -1;

	while(getline(infile, line))
	{
		if (line[0] == ';' || line[0] == '@')
		{
			continue;
		}
		if(line[0] == '>')
		{
			if (!seq.empty())
			{
				assert(copy_number != -1);
#ifdef ERROR
				if(seq.size() < K)
				{
					cerr << "[Info] Illegal input contigs, ignoring..." << endl;
					continue;
				}
#endif

				Edge edge(cur_index, seq.size(), seq, copy_number);

#ifdef TRAINING
				if (copy_number == 1)
				{
					unique_contigs_for_training.push(seq);
				}
#endif
				all_cdbg_edges.push_back(edge);
				++cur_index;
				seq.clear();
			}

			fields = Split_String(line).split("\t ");
			copy_number = atoi(fields[1].c_str());
		}
		else
		{
			seq += line;
		}
	}
	if (!seq.empty())
	{
		all_cdbg_edges.push_back(Edge( cur_index , seq.size() , seq , copy_number));
		seq.clear();
	}
	cout << "[Info] all contig number = " << all_cdbg_edges.size() << endl;
	infile.close();
	infile.clear();

#ifdef TRAINING
	//sort(uniq_edges_for_training.begin(), uniq_edges_for_training.end(),myFunc);

	cout << "[Info] Output unique contigs for training..." << endl;

	char buf[10];
	buf[0] = '\0';
	sprintf(buf, "%d",K);
	string fileName;
	fileName += string(buf);
	fileName += string("mer");
	fileName += ".unique_contig_for_training.fasta";

	ofstream out(fileName.c_str());
	for(int i = 0 ; i < MIN_EDGE_COUNT_FOR_TRAINING && !unique_contigs_for_training.empty(); ++i)
	{
		out << ">seq_" << i << endl;
		out << unique_contigs_for_training.top() << endl;
        unique_contigs_for_training.pop();
	}
    while(!unique_contigs_for_training.empty())
    {
        unique_contigs_for_training.pop();
    }
	out.close();
	out.clear();
//	unique_contigs_for_training.clear();

#endif

	infile.open(initial_edge_file_name.c_str());
	cur_index = 0;
	copy_number = 0;
	string dummy;
	int temp;

	while(!infile.eof())
	{

		infile >> temp;
		infile >> temp;
		infile >> seq;

		if(infile.eof())break;

#ifdef ERROR
		if(seq.size() < K)
		{
			cout << seq << endl;
			cout << K << endl;
			cerr << "[Info] Input contig illegal, ignoring..." << endl;
			getline(infile,dummy);
			getline(infile,dummy);
			if(infile.eof())
				break;
			seq.clear();
			continue;
		}
#endif
		Edge edge(cur_index, seq.size(), seq, copy_number);
		all_initial_edges.push_back(edge);

		cur_index++;
		seq.clear();

		getline(infile,dummy);

		getline(infile,dummy);
		if(infile.eof())
			break;
	}

	cout << "[Info] All condensed edges number = " << all_initial_edges.size() << endl;

	infile.close();
	infile.clear();

	cout << "[Info] Read contigs end..." << endl;
}

//Run BFS to fill current gap 
void Gap_Filling::BFS( const  int gap_index, const  int left_index, const  int right_index, int dis)
{

	vector<vector<int> >  final_sets; //only contain the middle sequence except two ends
	choice_type prev_sets;
	choice_type next_sets;

	dis += (K-1+3*var); //distance constraints
	dis += 30;

	int step = 0;

	vector<int> next = next_candidates_in_cdbg[left_index];

	assert(next.size() == ALPHABET-1);

	for(int i = 0 ; i < next.size() ; ++i)
	{
		int index = next[i];
		if( index == -1) continue;

		int len = all_cdbg_edges[index].get_len() - (K-1);

		if ( index == right_index)
		{
			final_sets.push_back(vector<int>(1,index));
		}
		prev_sets.push_back(make_pair( vector< int > ( 1 , index ) , len )) ;

	}
	++step;


	while(!prev_sets.empty())
	{
		++step;
		if(step > STEP)
		{
			cout << "[Info] BFS step is bigger than STEP..." << endl;
			break;
		}

		if(prev_sets.size() > 1000) 
		{
			break;
		}

		for(int i = 0 ; i < prev_sets.size(); ++i)
		{
			vector<int> path = prev_sets[i].first;
			int index = path[path.size()-1];
			int len = prev_sets[i].second;
			if ( len > dis ) continue;

			vector<int> next = next_candidates_in_cdbg[index];

			for(int i = 0 ; i < next.size() ; ++i)
			{
				if (next[i] == -1) continue;
				int node = next[i];
				if(node == right_index)
				{
					final_sets.push_back(path);
				}
				int len2 = len + all_cdbg_edges[node].get_len()-(K-1);
				path.push_back(node);
				next_sets.push_back(make_pair(path,len2));
				path.pop_back();
			}
		}
		prev_sets.swap(next_sets);
		next_sets.clear();
	}

	if(!final_sets.empty())
	{
		if (final_sets.size() == 1)
		{
			uniq_gap_info.push_back(make_pair(gap_index,final_sets[0]));
		}
		else
		{
			assert(final_sets.size() > 1);
			if (final_sets.size() > MAX_CHOICE)
			{
				fail_gap_info.push_back(gap_index);
			}
			else
			{
				assert(final_sets.size() <= MAX_CHOICE);
				multi_gap_info.push_back(make_pair(gap_index, final_sets));
			}
		}

		assert(gap_state.find(gap_index) == gap_state.end());
		gap_state.insert(gap_index);
		prev_sets.clear();
		next_sets.clear();
		final_sets.clear();
		return ;
	}

	//run BFS on the initial edges for fail gap
	prev_sets.clear();
	next_sets.clear();

	assert(final_sets.empty());
	final_sets.clear();

	step = 0;

	string start = all_cdbg_edges[left_index].get_suffix();
	string end = all_cdbg_edges[right_index].get_prefix();

	set<int> ends;

	candidates_iterator it = tails_initial_hash.find(end);

	if(it == tails_initial_hash.end())
	{
		fail_gap_info.push_back(gap_index);
		return;
	}
	assert(it->second.size() == ALPHABET-1);

	for(int i = 0 ; i < it->second.size() ; ++i)
	{
		if(it->second[i] == -1) continue;
		else ends.insert(it->second[i]);
	}
	if(ends.empty())
	{
		fail_gap_info.push_back(gap_index);
		return;
	}

	//find the first step indexs in initial
	it = heads_initial_hash.find(start);

	if(it == heads_initial_hash.end())
	{
		fail_gap_info.push_back(gap_index);
		prev_sets.clear();
		next_sets.clear();
		final_sets.clear();

		return;
	}
	assert(it->second.size() == ALPHABET-1);

	for(int i = 0 ; i < it->second.size() ; ++i)
	{
		int index = (it->second)[i];
		if( index == -1) continue;

		int len = all_initial_edges[index].get_len() - (K-1);

		if ( ends.find(index) != ends.end() )
		{
			final_sets.push_back(vector<int>(1,index));
		}
		prev_sets.push_back(make_pair( vector< int > ( 1 , index ) , len )) ;
	}
	++step;


	while(!prev_sets.empty())
	{
		++step;
		if(step > STEP) break;
		if(prev_sets.size() > 2000) 
		{
			break;
		}
		for(int i = 0 ; i < prev_sets.size(); ++i)
		{
			vector<int> path = prev_sets[i].first;
			int index = path[path.size()-1];
			int len = prev_sets[i].second;

			if(len > dis)
			{
				continue;
			}

			vector<int> next = next_candidates_in_initial[index];
			assert(next.size() == ALPHABET-1);

			for(int i = 0 ; i < ALPHABET-1; ++i)
			{
				int node = next[i];
				if(node == -1) continue;

				//if(node == right_index)
				if ( ends.find(node) != ends.end() )
				{
					final_sets.push_back(path);
				}
				int len2 = len + all_initial_edges[node].get_len()-(K-1);
				path.push_back(node);
				next_sets.push_back(make_pair(path,len2));
				path.pop_back();
			}
		}
		prev_sets.swap(next_sets);
		next_sets.clear();
	}

	if (!final_sets.empty())
	{
		if (final_sets.size() == 1)
		{
			uniq_gap_info.push_back(make_pair(gap_index,final_sets[0]));
		}	
		else
		{
			if (final_sets.size() > MAX_CHOICE)
			{
				fail_gap_info.push_back(gap_index);
			}
			else
			{
				assert(final_sets.size() <= MAX_CHOICE);
				multi_gap_info.push_back(make_pair(gap_index, final_sets));
			}
		}

		prev_sets.clear();
		next_sets.clear();
		final_sets.clear();
		return;
	}

	fail_gap_info.push_back(gap_index);
	prev_sets.clear();
	next_sets.clear();
	final_sets.clear();
	return;
}

//read scaffolds
void Gap_Filling::set_scaffolds()
{

	ifstream infile(scaffold_file_name.c_str());
	if (!infile)
	{
		cerr << "[Info] " << scaffold_file_name << " No Such file or directory..." << endl;
		exit(EXIT_FAILURE);
	}
	string line;

	int cur_scaff_index(0);
	int edge_number(0);

	vector< int > edges;
	vector< int > distances;
	vector<string > temp_split;

	try{

	    while(getline( infile, line ) )
	    { 
		    edges.clear();
			distances.clear();

			if(line[0] == ';' || line.empty())
			{
				continue;
			}
			if(line[0] == '>')
			{
				edges.clear();

				getline(infile, line);
				temp_split = Split_String(line).split("\t ");
				
				if(temp_split.size() == 0)
					continue;

				for(int i = 0 ;i < temp_split.size() ; ++i)
					edges.push_back( atoi(temp_split[i].c_str()));

				edge_number = edges.size();

				getline(infile, line);
				temp_split = Split_String(line).split("\t ");

				for(int i = 0 ;i < temp_split.size() ; ++i)
					distances.push_back( atoi(temp_split[i].c_str()));
				assert( distances.size()+1 == edges.size());

				Scaffold scaff(cur_scaff_index, edge_number, edges, distances);
				all_scaffolds.push_back(scaff);

				++cur_scaff_index;
			}
		
		}
	}
    //Illegal input format detected
	catch(runtime_error err)
	{
		cout << err.what() << "Input Scaffolds Error!!!" << endl;
		infile.close();
		infile.clear();
		exit(1);
	}

	infile.close();
	infile.clear();
	cout << "[Info] Input Scaffolds Number = " << all_scaffolds.size() << endl;
}


//get gap info and classify gaps
void Gap_Filling::get_gap_info()
{
	cout << "[Info] Begin Gap Filling ... " << endl;

	set_edges();
	set_scaffolds();

	cout << "[Info] Set edges and Scaffolds end..." << endl;
	cout << "[Info] Begin build condensed de Buijn graph... " << endl;

	string prefix, suffix;
	candidates_iterator it;

	//hash initial node
	for(int i = 0 ; i < all_initial_edges.size(); ++i)
	{

		prefix = all_initial_edges[i].get_prefix();
		char next_base = all_initial_edges[i].get_pre_base();
		int index = -1;
		switch(next_base)
		{
			case 'A': index = A; break;
			case 'T': index = T; break;
			case 'C': index = C; break;
			case 'G': index = G; break;
					  //default G
			default : index = G; break;
		}
		assert(index != -1);

		it = heads_initial_hash.find(prefix);

		if ( it == heads_initial_hash.end())
		{
			heads_initial_hash.insert(make_pair(prefix,vector<int>(ALPHABET-1,-1)));
			heads_initial_hash[prefix][index] = i;
		}
		else
		{
			/*
			   if(heads_initial_hash[prefix][index] != -1)
			   {
			   cout << heads_initial_hash[prefix][index] << " " << i <<  "share the same head..." << endl;
			   }
			   */
			heads_initial_hash[prefix][index] = i;
		}

		suffix = all_initial_edges[i].get_suffix();
		next_base = all_initial_edges[i].get_suf_base();
		index = -1;
		switch(next_base)
		{
			case 'A': index = A; break;
			case 'T': index = T; break;
			case 'C': index = C; break;
			case 'G': index = G; break;
					  //default G
			default : index = G; break;
		}
		assert(index != -1);

		it = tails_initial_hash.find(suffix);

		if (it == tails_initial_hash.end())
		{

			tails_initial_hash.insert(make_pair(suffix,vector<int>(ALPHABET-1,-1)));
			tails_initial_hash[suffix][index] = i;
		}
		else
		{
			/*
			   if(tails_initial_hash[suffix][index] != -1)
			   {
			   cout << tails_initial_hash[suffix][index] << " " << i << "share the same tail..." << endl;
			   }
			   */
			tails_initial_hash[suffix][index] = i;
		}
	}

	cout << "[Info] Begin build the prev and next node of contig in CDBG..." << endl;

	prev_candidates_in_initial.resize(all_initial_edges.size());
	next_candidates_in_initial.resize(all_initial_edges.size());

	for(int i = 0 ; i < all_initial_edges.size();++i)
	{
		suffix = all_initial_edges[i].get_suffix();

		it = heads_initial_hash.find(suffix);
		next_candidates_in_initial[i]= vector<int>(ALPHABET-1,-1);

		if(it != heads_initial_hash.end())
		{
			assert(it->second.size() == ALPHABET-1);
			for(int k = 0 ; k < ALPHABET-1 ; ++k)
			{
				next_candidates_in_initial[i][k] = (it->second[k]);
			}
		}

		string prefix = all_initial_edges[i].get_prefix();
		prev_candidates_in_initial[i] = vector<int>(ALPHABET-1,-1);

		it = tails_initial_hash.find(prefix);
		if (it != tails_initial_hash.end())
		{
			assert(it->second.size() == ALPHABET-1);
			for(int k = 0 ; k < ALPHABET-1; ++k)
			{
				prev_candidates_in_initial[i][k] = (it->second[k]);
			}
		}
	}

	cout << "[Info] Begin build condensed de Buijn graph... " << endl;

	//hash cdbg node
	for(int i = 0 ; i < all_cdbg_edges.size(); ++i)
	{
		prefix = all_cdbg_edges[i].get_prefix();
		char next_base = all_cdbg_edges[i].get_pre_base();
		int index = -1;
		switch(next_base)
		{
			case 'A': index = A; break;
			case 'T': index = T; break;
			case 'C': index = C; break;
			case 'G': index = G; break;
					  //default G
			default : index = G; break;
		}

		it = heads_cdbg_hash.find(prefix);
		if ( it == heads_cdbg_hash.end())
		{
			heads_cdbg_hash.insert(make_pair(prefix,vector<int>(ALPHABET-1,-1)));
			heads_cdbg_hash[prefix][index] = i;
		}
		else
		{
			//assert(heads_cdbg_hash[prefix][index] == -1);
			heads_cdbg_hash[prefix][index] = i;
		}

		suffix = all_cdbg_edges[i].get_suffix();
		next_base = all_cdbg_edges[i].get_suf_base();
		index = -1;
		switch(next_base)
		{
			case 'A': index = A; break;
			case 'T': index = T; break;
			case 'C': index = C; break;
			case 'G': index = G; break;
					  //default G
			default : index = G; break;
		}

		it = tails_cdbg_hash.find(suffix);
		if ( it == tails_cdbg_hash.end())
		{

			tails_cdbg_hash.insert(make_pair(suffix,vector<int>(ALPHABET-1,-1)));
			tails_cdbg_hash[suffix][index] = i;
		}
		else
		{
			//assert(tails_cdbg_hash[suffix][index] == -1);
			tails_cdbg_hash[suffix][index] = i;
		}
	}

	cout << "[Info] Set edges in condensed de Bruijn graph..." << endl;
	next_candidates_in_cdbg.resize(all_cdbg_edges.size());
	prev_candidates_in_cdbg.resize(all_cdbg_edges.size());

	for(int i = 0 ; i < all_cdbg_edges.size();++i)
	{
		next_candidates_in_cdbg[i] = vector<int>(ALPHABET-1,-1);

		suffix = all_cdbg_edges[i].get_suffix();
		it = heads_cdbg_hash.find(suffix);

		if(it != heads_cdbg_hash.end())
		{
			assert(it->second.size() == ALPHABET-1);
			for(int k = 0 ; k < ALPHABET-1 ; ++k)
			{
				next_candidates_in_cdbg[i][k] = (it->second[k]);
			}
		}

		prev_candidates_in_cdbg[i] = vector<int>(ALPHABET-1,-1);

		string prefix = all_cdbg_edges[i].get_prefix();
		it = tails_cdbg_hash.find(prefix);

		if (it != tails_cdbg_hash.end())
		{
			assert(it->second.size() == ALPHABET-1);
			for(int k = 0 ; k < ALPHABET-1; ++k)
			{
				prev_candidates_in_cdbg[i][k] = (it->second[k]);
			}
		}
	}

	cout << "[Info] Build condensed de Bruijn graph end..." << endl;

	assert(next_candidates_in_initial.size() == all_initial_edges.size());
	assert(next_candidates_in_cdbg.size() == all_cdbg_edges.size());
	assert(prev_candidates_in_initial.size() == all_initial_edges.size());
	assert(prev_candidates_in_cdbg.size() == all_cdbg_edges.size());

	heads_cdbg_hash.clear();
	tails_cdbg_hash.clear();


	set< pair<int, int> > visited; //record the processed gaps
	int cur_gap_index(0);

	for(vector<Scaffold >::size_type i = 0 ; i < all_scaffolds.size(); ++i)
	{
		if (all_scaffolds[i].get_edge_numbers() == 1)
		{
			continue;
		}
		vector<int > edge_indexs = all_scaffolds[i].get_edges_of_scaffold();

		for(int j = 0 ; j <  edge_indexs.size()-1; ++j )
		{
			int left_index = edge_indexs.at(j);
			int right_index = edge_indexs.at(j+1);

			gap_indexs.push_back(make_pair(left_index,right_index));
			gap_distances.push_back(all_scaffolds[i].get_specific_pos_distance_of_scaffold(j));

			string suffix = all_cdbg_edges.at(left_index).get_suffix();
			string prefix = all_cdbg_edges.at(right_index).get_prefix();

			int temp_overlap = 0;
			if( ( temp_overlap = alignment(suffix, prefix)) >= overlap)
			{
				assert(visited.find(make_pair(left_index, right_index)) == visited.end());
				visited.insert(make_pair(left_index, right_index));
				pre_gap_info.insert(make_pair(cur_gap_index,temp_overlap));
			}
			cur_gap_index++;
		}
	}
	cout << "[Info] The number of gaps = " << gap_indexs.size() << endl;
	cout << "[Info] Filling the gaps which shares overlap ( < " << overlap << ") = " << pre_gap_info.size() << endl;

	assert(gap_indexs.size() == cur_gap_index);
	assert(gap_indexs.size() == gap_distances.size());

#ifdef RESULT 
    int count_base_in_gap=0;
#endif

	for(int i = 0 ; i < gap_indexs.size(); ++i)
	{
		gap_indexs_hash.insert(make_pair(gap_indexs[i],i));
	}

	assert(gap_indexs.size() == gap_indexs_hash.size());

	//for each gap , run BFS and get uniq_gap_info and multi_gap_info

	for(int i = 0 ; i < gap_indexs.size(); ++i )
	{

		int left_index = gap_indexs[i].first;
		int right_index = gap_indexs[i].second;

		if (visited.find(gap_indexs[i]) != visited.end())
		{
			continue;
		}
		int distance = gap_distances[i];

		if (distance > ( MU + 3*var ) )
		{
			fail_gap_info.push_back(i);
			continue;
		}
		try
		{
			BFS(i, left_index, right_index, distance);
#ifdef RESULT 
            count_base_in_gap += distance;
#endif
            
		}
		catch(bad_alloc)
		{
			cerr << "[Info] BFS is too memory-intensive, ignoring..." << endl;
			fail_gap_info.push_back(i);
		}
	}

#ifdef  RESULT
    cout << "[Info] Base counts in gaps : " << count_base_in_gap << endl;
#endif

	cout << "[Info] Unique candidate gaps number = " << uniq_gap_info.size() << endl;
	cout << "[Info] Multi candidates gaps number = " << multi_gap_info.size() << endl;
	cout << "[Info] Failed gaps number = " << fail_gap_info.size() << endl;
	next_candidates_in_cdbg.clear();
	prev_candidates_in_cdbg.clear();
	next_candidates_in_initial.clear();
	prev_candidates_in_initial.clear();
	heads_initial_hash.clear();
	tails_initial_hash.clear();
}

//Filling the unique candidate gaps
//contain two ends
void Gap_Filling::uniq_candidate_gap_filling()
{
	cout << "[Info] Start filling unique candidates gaps ... " << endl;
	int uniq_gap_number = uniq_gap_info.size();

	string seq;

	for(int i = 0 ;  i < uniq_gap_number ;  ++i )
	{
		int cur_gap_index = uniq_gap_info[i].first;
		vector<int> &path = uniq_gap_info[i].second;

		int left_index = gap_indexs[cur_gap_index].first;
		int right_index = gap_indexs[cur_gap_index].second;

		seq.clear();

		//filled by condensed contigs
		if (gap_state.find(cur_gap_index) != gap_state.end())
		{
            seq.clear();
			seq += all_cdbg_edges[left_index].get_edge_seq();
			assert(seq.size() >= K);

			seq.resize(seq.size()-(K-1));

			for( int j = 0; j < path.size(); ++j)
			{
				assert(all_cdbg_edges[path[j]].get_len() >= K);
				seq += (all_cdbg_edges[path[j]].get_edge_seq());
				assert(seq.size() >= K);
				seq.resize(seq.size()-(K-1));
			}
			seq += (all_cdbg_edges[right_index].get_edge_seq());
			uniq_gap_seq.insert(make_pair(cur_gap_index,seq));
		}
		//filled by initial contigs
		else
		{
            seq.clear();
			assert(gap_state.find(cur_gap_index) == gap_state.end());

			seq += all_cdbg_edges[left_index].get_edge_seq();
			assert(seq.size() >= K);
			seq.resize(seq.size()-(K-1));

			for( int j = 0; j < path.size(); ++j)
			{
				assert(all_initial_edges[path[j]].get_len() >= K);
				seq += (all_initial_edges[path[j]].get_edge_seq());
				assert(seq.size() >= K);
				seq.resize(seq.size()-(K-1));
			}
			seq += (all_cdbg_edges[right_index].get_edge_seq());
			uniq_gap_seq.insert(make_pair(cur_gap_index,seq));
		}
	}
	assert(uniq_gap_seq.size() == uniq_gap_info.size());
	cout << "[Info] Fill unique candidate gaps end..." << endl;
}

//Do not contain the two ends
void Gap_Filling::multi_candidates_gap_filling()
{
	cout << "[Info] Start multi candidates gap filling... " << endl;

	string seq;
	vector<string> seq_s;

	for(int i = 0; i < multi_gap_info.size(); ++i)
	{
		assert(multi_gap_info[i].second.size() <= MAX_CHOICE);

		int cur_gap_index = multi_gap_info[i].first;
		int left_index = gap_indexs[cur_gap_index].first;
		int right_index = gap_indexs[cur_gap_index].second;

		seq_s.clear();

		vector<vector<int> > & paths = multi_gap_info[i].second;

		assert(paths.size() > 1);

		//The gap is filled by contigs in cdbg 
		if(gap_state.find(cur_gap_index) != gap_state.end())
		{
			seq.clear();
            seq_s.clear();
			for( int j = 0 ; j < paths.size(); ++j)
			{
				seq.clear();
				seq += all_cdbg_edges[left_index].get_edge_seq();
                assert(seq.size() >= K);
			    seq.resize(seq.size()-(K-1));

				for( int k = 0; k < paths[j].size(); ++k)
				{
					seq += (all_cdbg_edges[paths[j][k]].get_edge_seq());
					assert(seq.size() >= K);
					seq.resize(seq.size()-(K-1));
				}
				seq += (all_cdbg_edges[right_index].get_edge_seq());
				seq_s.push_back(seq);
			}
			multi_gap_seq.insert(make_pair(cur_gap_index, seq_s));
		}
		else
		{
			//The gap is filled by contigs in initial CDBG
			seq.clear();
            seq_s.clear();
			for( int j = 0 ; j < paths.size(); ++j)
			{
				seq.clear();
				seq += all_cdbg_edges[left_index].get_edge_seq();
                assert(seq.size() >= K);
                seq.resize(seq.size()-(K-1));

				for( int k = 0; k < paths[j].size(); ++k)
				{
					assert(paths[j][k] < all_initial_edges.size());

					seq += (all_initial_edges[paths[j][k]].get_edge_seq());
					assert((all_initial_edges[paths[j][k]].get_edge_seq()).size() >= K);
					seq.resize(seq.size()-(K-1));		
				}
				seq += (all_cdbg_edges[right_index].get_edge_seq());
				seq_s.push_back(seq);
			}
			multi_gap_seq.insert(make_pair(cur_gap_index, seq_s));
		}
	} 

	assert(multi_gap_seq.size() == multi_gap_info.size());

	cout << "[Info] Multi candidate gap filling end ... " << endl;

	//clear the initial edges since it would not be used again
	all_initial_edges.clear();
}

/*
 * build hash table  map<pair<int, int> , int >
 * uniq -- 1
 * multi --  >1
 * pre -- -(overlap)
 * fail -- 0
 * */

void Gap_Filling::get_gap_candidate_number_hash()
{

	//unique gap -- 1
	for(int i = 0 ; i < uniq_gap_info.size() ; ++i)
	{
		pair<int,int> edge_pair = gap_indexs[uniq_gap_info[i].first];
		gap_candidate_number_hash.insert(make_pair(edge_pair, 1));
	}
	assert(gap_candidate_number_hash.size() == uniq_gap_info.size());

	//multi gap -- #choices
	for(int i = 0 ; i < multi_gap_info.size(); ++i)
	{
		pair<int,int> edge_pair = gap_indexs[multi_gap_info[i].first];
		int candidate_number = multi_gap_info[i].second.size();

		assert(candidate_number <= MAX_CHOICE);
		gap_candidate_number_hash.insert(make_pair(edge_pair,candidate_number));
	}
	assert(gap_candidate_number_hash.size() == uniq_gap_info.size()+ multi_gap_info.size());

	//pre gap info -- -(overlap)
	for(map<int,int>::iterator it = pre_gap_info.begin(); it != pre_gap_info.end(); ++it)
	{
		pair<int,int> edge_pair = gap_indexs[it->first];
		gap_candidate_number_hash.insert(make_pair(edge_pair,-(it->second)));
	}
	//fail gap info -- 0
	for(int i = 0; i < fail_gap_info.size() ;++i)
	{
		pair<int,int> edge_pair = gap_indexs[fail_gap_info[i]];
		gap_candidate_number_hash.insert(make_pair(edge_pair,0));
	}

	assert(gap_candidate_number_hash.size() == uniq_gap_info.size()+ multi_gap_info.size() + pre_gap_info.size() + fail_gap_info.size());
	assert(gap_candidate_number_hash.size() == gap_indexs.size());
}

bool myFuncPairIntLong(pair<int,long> lhs, pair<int,long> rhs)
{
	return lhs.second > rhs.second;
}

//get initial scaffolds with no filling gap filling and using 'N' to fill failed gaps
void Gap_Filling::output_initial_scaffolds_seq()
{
	//get gap number hash
	get_gap_candidate_number_hash();
	int cur_gap_index = 0 ;

	long GENOME = 0;		//estimate the genome size using total length of scaffold

	char buf[10];
	buf[0] = '\0';
	sprintf(buf, "%d",K);

	string fileName;
	fileName += string(buf);
	fileName += string("mer");
	fileName += ".scaf_seq_with_gaps";

	//(index , length) pair of scaffold
	vector<pair<int, long> > scaf_index_len ;
	ofstream out_scaf_seq(fileName.c_str());
	if(!out_scaf_seq)
	{
		cerr << "[Info] Create " << fileName << " failed..." << endl;
		exit(EXIT_FAILURE);
	}

	string seq;         //scaffold sequence
	string one_segment; //one segment of scaffold since the multi gaps

	//for scaffolds segments
	vector<string> cur_scaff_segments;

	map<pair<int,int>, int>::iterator it;
	for(int i = 0 ; i < all_scaffolds.size(); ++i)
	{
		seq.clear();
		one_segment.clear();
		cur_scaff_segments.clear();

		vector<int> edge_indexs = all_scaffolds[i].get_edges_of_scaffold();
		assert(edge_indexs.size() > 0);
		seq += all_cdbg_edges[edge_indexs[0]].get_edge_seq();
		one_segment += all_cdbg_edges[edge_indexs[0]].get_edge_seq();

		if(edge_indexs.size() == 1)
		{
			scaf_index_len.push_back(make_pair(i,seq.size()));
			out_scaf_seq << ">scaf_" << i << "_" << seq.size()  << endl;
			out_scaf_seq << seq << endl;
			GENOME += seq.size();

			//for scaf segment
			cur_scaff_segments.push_back(one_segment);
			all_scaffolds_segments.push_back(cur_scaff_segments);
			continue;
		}

		for(int j=0 ; j+1 < edge_indexs.size(); ++j)
		{
			int left_index = edge_indexs[j];
			int right_index = edge_indexs[j+1];

			int left_len = all_cdbg_edges[left_index].get_len();
			int right_len = all_cdbg_edges[right_index].get_len();

			pair<int,int> edge_pair = make_pair(left_index,right_index);

			assert(gap_indexs[cur_gap_index].first == left_index && gap_indexs[cur_gap_index].second == right_index);

			it	= gap_candidate_number_hash.find(edge_pair);

			assert( it != gap_candidate_number_hash.end());
			//overlap
			if(it->second < 0)
			{
				assert(one_segment.size() > abs(it->second));
				one_segment.resize(one_segment.size() + it->second);
				one_segment += all_cdbg_edges[right_index].get_edge_seq();
			}
			//unique
			else if(it->second == 1 )
			{
				assert(one_segment.size() >= left_len);
				one_segment.resize(one_segment.size() - left_len);
				one_segment += uniq_gap_seq[cur_gap_index];
			}
			//failed gap
			else if(it->second == 0)
			{
				//for failed gap , the count of 'N' is determined by the distance estimation
				if(gap_distances[cur_gap_index] > (MU+3*var))
				{
					one_segment += string(MU+3*var,'N');
				}
				else if (gap_distances[cur_gap_index] < 0)
				{
					;
				}
				else
				{
					one_segment += string(gap_distances[cur_gap_index],'N');
				}
				one_segment += all_cdbg_edges[right_index].get_edge_seq();
			}
			//multi_gap
			else
			{
				assert(it->second <= MAX_CHOICE);
				cur_scaff_segments.push_back(one_segment);
				initial_multi_gaps.push_back(cur_gap_index);
				seq += one_segment;
				one_segment.clear();
				one_segment += all_cdbg_edges[right_index].get_edge_seq();
			}
			cur_gap_index++;
		}
		cur_scaff_segments.push_back(one_segment);
        assert(cur_scaff_segments.size() > 0 );

		all_scaffolds_segments.push_back(cur_scaff_segments);

		//initial scaff
		seq += one_segment;
		out_scaf_seq << ">scaf_" << i << "_" <<  seq.size()  << endl;
		out_scaf_seq << seq << endl;
		GENOME += seq.size();
		scaf_index_len.push_back(make_pair(i,seq.size()));
	}
	assert(cur_gap_index == gap_candidate_number_hash.size());
	assert(all_scaffolds_segments.size() ==  all_scaffolds.size());
	assert(all_scaffolds.size() ==  scaf_index_len.size());
	assert(initial_multi_gaps.size() == multi_gap_info.size());

	out_scaf_seq.close();
	out_scaf_seq.clear();


	//output result
	fileName.clear();
	fileName += string(buf);
	fileName += string("mer");
	fileName += ".scaf_len";

	ofstream out_scaf_len(fileName.c_str());
	if(!out_scaf_len)
	{
		cerr << "[Info] " << fileName << " No such file or directory!!!" << endl;
		exit(EXIT_FAILURE);
	}

    fileName.clear();
    fileName += string(buf);
    fileName += string("mer");
    fileName += ".final_result";

	ofstream out_final_result(fileName.c_str());
	if(!out_final_result)
	{
		cerr << "[Info] " << fileName << " No such file or directory!!!" << endl;
		exit(EXIT_FAILURE);
	}

	//find top 20
	sort(scaf_index_len.begin(), scaf_index_len.end(), myFuncPairIntLong);

    for(int i=0;i < scaf_index_len.size();++i)
    {
        out_scaf_len << scaf_index_len[i].second << endl;
    }
    out_scaf_len.close();
    out_scaf_len.clear();

	out_final_result << "Genome estimated size\t:\t" << GENOME << endl;

	long N50 = GENOME * 0.5 ;
	long N90 = GENOME * 0.9 ;

	out_final_result << "---Scaffold---"<< endl;
	out_final_result << "Scaff length : Top 20..." << endl;

	int size = scaf_index_len.size();
	for(int i = 0 ; i < min(20, size); ++i)
	{
		out_final_result << scaf_index_len[i].second << endl;
	}

	vector<long> sum (scaf_index_len.size(),0);
	
	sum[0] = scaf_index_len[0].second;
	for(int i = 1; i < scaf_index_len.size(); ++i)
	{
		sum[i] +=  sum[i-1] + scaf_index_len[i].second;
	}
	long N50_val = 0;
	long N90_val = 0;

	for(int i = 0 ; i < sum.size(); ++i)
	{
		if(sum[i] >= N50 && N50_val == 0)
		{
			N50_val = scaf_index_len[i].second;
		}
		if (sum[i] >= N90 && N90_val == 0)
		{
			N90_val = scaf_index_len[i].second;
			break;
		}
	}

	out_final_result << "N50\t:\t" << N50_val << endl;
	out_final_result << "N90\t:\t" << N90_val << endl;

	out_final_result << "All gaps number in scaffolds\t:\t" << gap_indexs.size() << endl;
	out_final_result << "filled gap numbers\t:\t" << gap_indexs.size() - fail_gap_info.size() << endl;
	out_final_result << "failed gap numbers\t:\t" << fail_gap_info.size() << endl;

	out_final_result.close();
	out_final_result.clear();

}

//DFS to get the indexs of new multi gaps to vetify the corresponding index of initial gaps
void Gap_Filling::DFS(
		vector<int> & edge_indexs, 
		const vector<int> &merged, 
		vector<string> &indexs, 
		int pos , 
		vector<string> &ret_seqs, 
		string index_path, 
		string &path , 
		vector<string> &adjacent)
{
	//the end
	assert(adjacent[pos].size() <= MU+var*3);
	path += adjacent[pos];
	if (pos == merged.size())
	{
		ret_seqs.push_back(path);
		index_path.resize(index_path.size()-1);
		indexs.push_back(index_path);
		return;
	}
	string temp = index_path;
	string temp_seq = path;
	for(int i = 0 ; i < merged[pos]; ++i)
	{
		index_path += ('0' + i);
		index_path += "_";

		path.resize(path.size() - (K-1));
		path += multi_gap_seq[edge_indexs[pos]][i];

		DFS(edge_indexs, merged, indexs, pos+1, ret_seqs, index_path, path , adjacent);
		index_path = temp;
		path = temp_seq;
	}

}


//For multi candidate gaps , produce local sequences to scoring
void Gap_Filling::get_multi_candidates_gap_local_seq()
{
    //DEBUG
    int debug_count = 0;
    for(int i=0;i<all_scaffolds_segments.size() ;++i)
    {
        debug_count += (all_scaffolds_segments[i].size()-1);
    }
    assert(debug_count == multi_gap_info.size());

	//clear useless data structure
	uniq_gap_info.clear();
	uniq_gap_seq.clear();

    //cout << all_scaffolds_segments.size() << endl;

	cout << "[Info] Output multi candidates gap local seq to scoring..." << endl;

	char buf[10];
	buf[0] = '\0';
	sprintf(buf, "%d",K);

	string fileName;
	fileName += string(buf);
	fileName += "mer";
	fileName += ".multi_candidates_seq_to_PerM.fasta";

	ofstream out_multi_gaps(fileName.c_str());
	if (!out_multi_gaps)
	{
		cerr << "[Info] Create " << fileName << " No such file or directory!!!" << endl;
		exit(EXIT_FAILURE);
	}

	buf[0] = '\0';
	sprintf(buf, "%d",K);
	fileName.clear();
	fileName += string(buf);
	fileName += "mer";
	fileName += ".multi_gap_count_to_score";

	//ofstream outfile2((TEMP_FILE_DIR + fileName).c_str());
	ofstream out_gap_count(fileName.c_str());
	if (!out_gap_count)
	{
		cerr << "[Info] " << fileName << " No such file or directory!!!" << endl;
		exit(EXIT_FAILURE);
	}

	buf[0] = '\0';
	sprintf(buf, "%d",K);
	fileName.clear();
	fileName += string(buf);
	fileName += "mer";
	fileName += ".multi_scaff_segment_to_score";
	ofstream out_scaff_segment(fileName.c_str());
	if(!out_scaff_segment)
	{
		cerr << "[Info] " << fileName << " No such file or directory!!!" << endl;
		exit(EXIT_FAILURE);
	}


	srand(time(NULL));

	int threshold = MU + var*3;
	int cur_gap_index = 0;				
	int index = 0;				//record the current multi gap index
	int output_gap_index = 0;

	vector<string> candidate_seqs;		//store current gap's all local sequences

	vector<int> merged_gap_choices;		//The gap candidates need to be merged
	vector<int> merged_gap_indexs;		//The gaps need to be merged

	vector<string> cur_scaffold_segment;

	vector<string> adjacent;

	for(int i = 0 ; i < all_scaffolds_segments.size() ; ++i)
	{
		cur_scaffold_segment.clear();

		vector<string> seqs = all_scaffolds_segments[i];

		assert(seqs.size() >= 1);
		if(seqs.size() == 1)
		{
			out_scaff_segment << ">seq\t" << i << endl;
			out_scaff_segment << seqs[0] << endl; 
			continue;
		}
		int j = 1, cur = 0;
		int merge_gap_count = 0;

        cur_scaffold_segment.push_back(seqs[0]);
		while(j < seqs.size())
		{
			merge_gap_count = 0;
			while(j < seqs.size()-1 && merge_gap_count <= MAX_MERGE_GAP && seqs[j].size() < threshold)
			{
				++j;
				++merge_gap_count;
			}

			merged_gap_choices.clear();
			merged_gap_indexs.clear();
			candidate_seqs.clear();
			adjacent.clear();

			for(int k = cur; k < j; ++k)
			{
				merged_gap_indexs.push_back(initial_multi_gaps[index]);
				merged_gap_choices.push_back((multi_gap_seq[initial_multi_gaps[index]]).size());
				adjacent.push_back(seqs[k]);
				++index;
			}

			adjacent.push_back(seqs[j]);

			assert(adjacent.size() == merged_gap_indexs.size()+1);
			assert(merged_gap_indexs.size() == (j-cur));
			assert(merged_gap_choices.size() == (j-cur));
/*
            if(merge_gap_count)
            {
              cout << "merge " << merge_gap_count << " seqs"  << endl;
              for(int z=0;z<adjacent.size();++z)
              {
                cout << "merge " << adjacent[z].size() << endl;
              }
              cout << "current index = " << index  << endl;
              cout << endl;
            }
*/

            //if current gaps candidates is greater MAX_CHOICE, then discard these multi-gaps and random choose the gap sequence
            long choices = 1;
            string temp_segment ;
            for(int i_count = 0 ; i_count < merged_gap_choices.size(); ++i_count)
            {
                choices *= merged_gap_choices[i_count];
            }

            if(choices > MAX_CHOICE)
            {
                int i_count = 0 ;
                for(; i_count < merged_gap_indexs.size(); ++i_count)
                {
                    //int cur_index = initial_multi_gaps[merged_gap_indexs[i_count]];
                    int cur_index = merged_gap_indexs[i_count];
                    temp_segment += adjacent[i_count];
                    assert(multi_gap_seq.find(cur_index)!= multi_gap_seq.end());
                    int rand_number = rand()%(multi_gap_seq[cur_index].size());
                    temp_segment += multi_gap_seq[cur_index][rand_number];
                    //temp_segment += multi_gap_seq[initial_multi_gaps[cur_index]][rand()%(multi_gap_seq[initial_multi_gaps[cur_index]].size())];
                }
                temp_segment += adjacent[i_count];
              
                cur_scaffold_segment.push_back(temp_segment);
    			cur = j;
    			j++;
                continue;
            }

			//Record the index sequence,as 0_0_1
			vector<string> indexs;
			string index_path;
			string path;

			//extract the prefix and suffix sequence of MU+3*var
			if (adjacent[0].size() > threshold)
			{
				adjacent[0] = adjacent[0].substr(adjacent[0].size()-threshold);
			}
			if (adjacent[adjacent.size()-1].size() > threshold)
			{
				adjacent[adjacent.size()-1].resize(threshold);
			}

			DFS(merged_gap_indexs, merged_gap_choices, indexs, 0, candidate_seqs, index_path, path, adjacent);

			int left_index = gap_indexs[merged_gap_indexs[0]].first;
			int right_index = gap_indexs[merged_gap_indexs[merged_gap_indexs.size()-1]].second;
			new_multi_gaps.push_back(make_pair(left_index, right_index));

			for(int l = 0 ; l < merged_gap_indexs.size(); ++l)
			{
				assert(multi_gap_candidate_seq.find(merged_gap_indexs[l]) == multi_gap_candidate_seq.end());
				multi_gap_candidate_seq[merged_gap_indexs[l]] = candidate_seqs;
			}
			//new_multi_gap_candidate_seq.insert(make_pair(new_gap_indexs.size()-1, candidate_seqs));
			out_gap_count << candidate_seqs.size() << endl;
			for(int m = 0 ; m < candidate_seqs.size(); ++m)
			{
				out_multi_gaps << ">seq_" << output_gap_index << "_" << m << endl;
				out_multi_gaps << candidate_seqs[m] << endl;
			}
			output_gap_index++;
			candidate_seqs.clear();
			cur = j;
			j++;
			cur_scaffold_segment.push_back(seqs[cur]);
		}
		out_scaff_segment << ">seq  " << i << endl;
		for(int j=0; j < cur_scaffold_segment.size(); ++j)
		{
			out_scaff_segment << cur_scaffold_segment[j] << endl;
		}
	}
}


//operator<<
ostream &operator<<(ostream &os, const Gap_Filling &obj )
{
	os << "Pre_gap_info..." << endl;
	for( map<int,int>::const_iterator it = (obj.pre_gap_info).begin(); it != (obj.pre_gap_info).end(); ++it)
	{
		os << "("
			<< (obj.gap_indexs)[it->first].first
			<< ","
			<<(obj.gap_indexs)[it->first].second
			<< ")"
			<< it->second 
			<< endl;
	}

	os << "uniq_gap_info..." << endl;
	for(int  i = 0 ; i < obj.uniq_gap_info.size() ; ++i)
	{	
		os << obj.gap_indexs.at(obj.uniq_gap_info[i].first).first
			<< " "
			<< obj.gap_indexs.at(obj.uniq_gap_info[i].first).second
			<< " : ";
		copy((obj.uniq_gap_info[i].second).begin(), (obj.uniq_gap_info[i].second).end(),
				ostream_iterator<int>(os,"\t"));
		os << endl;
	}

	os << "multi_gap_info..." << endl;
	for(int i = 0 ; i < obj.multi_gap_info.size(); ++i)
	{
		os << "( " 
			<< obj.gap_indexs.at(obj.multi_gap_info[i].first).first
			<< ","
			<< obj.gap_indexs.at(obj.multi_gap_info[i].first).second
			<< " ) "
			<< " : ";
		for(int j = 0 ; j < (obj.multi_gap_info[i].second).size(); ++j)
		{
			os << j << " = " ;
			copy( ((obj.multi_gap_info[i]).second)[j].begin(), ((obj.multi_gap_info[i]).second)[j].end(),ostream_iterator<int> (os,"\t"));
			os << endl;
		}
	}

	os << "uniq_gap_seq\n " << endl;
	for(map<int,string>::const_iterator it = obj.uniq_gap_seq.begin();
			it != obj.uniq_gap_seq.end();
			++it)
	{
		os << "( "
			<< obj.gap_indexs.at(it->first).first
			<< ","
			<< obj.gap_indexs.at(it->first).second
			<< " )\n"
			<< it->second
			<< endl;
	}

	os << "multi_gap_seq\n" << endl;
	int i = 0;

	for(map<int,vector<string> >::const_iterator it = obj.multi_gap_seq.begin();
			it != obj.multi_gap_seq.end();
			++it)
	{

		os << "( "
			<< obj.gap_indexs.at(it->first).first
			<< ","
			<< obj.gap_indexs.at(it->first).second
			<< " )\n";
		copy( (it->second).begin(), (it->second).end(), ostream_iterator<string>(os, "\n"));
		os << endl;
	}
	os << "each scaffold segment..." << endl;

	int sum = 0 ;
	int index = 0;
	assert(obj.all_scaffolds_segments.size() == obj.all_scaffolds.size());
	for(int i = 0; i < obj.all_scaffolds_segments.size() ; ++i)
	{
		assert(obj.all_scaffolds_segments[i].size() > 0);

		if (obj.all_scaffolds_segments[i].size() > 1)
		{
			sum += (obj.all_scaffolds_segments[i].size()-1);
		}
		os << "i = " << i << endl;

		int j = 0;
		for(; j < obj.all_scaffolds_segments[i].size()-1; ++j)
		{
		} index += j;

		copy(obj.all_scaffolds_segments[i].begin(), obj.all_scaffolds_segments[i].end(),ostream_iterator<string>(os,"\n"));
		os << endl;
	}

	index = 0;
	for( vector<vector<string> >::const_iterator 
			it = (obj.all_scaffolds_segments).begin();
			it != (obj.all_scaffolds_segments).end();
			++it)
	{

		copy((*it).begin(), (*it).end(),  ostream_iterator<string>(os,"\n"));
		for(int i = 0;i < (*it).size()-1; ++i,++index)
		{
			os 
				<< "current gap = " 
				<< (obj.gap_indexs)[i].first
				<< " " 
				<< (obj.gap_indexs)[i].second 
				<< endl;
		}
	}

	cout << "[Info] Output each scaffold segments end..." << endl;
	return os;
}
