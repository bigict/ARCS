
#include "Gap_Filling.h"

#define TRAINING
//#define DEBUG
//#define ERROR

Gap_Filling::Gap_Filling(const string & _contig_file_name, const string &_scaffold_file_name, const string  &_initial_contig_file_name)
	:edge_file_name(_contig_file_name),scaffold_file_name(_scaffold_file_name),initial_edge_file_name(_initial_contig_file_name)
{
	ifstream in(_contig_file_name.c_str());
	if (!in) 
	{
		cerr << _contig_file_name << " : No such file or directory!!!" << endl; 
		exit(1);
	}
	in.close();
	in.clear();

	in.open(_scaffold_file_name.c_str());
	if (!in) 
	{
		cerr << _scaffold_file_name << " : No such file or directory!!!" << endl;
		exit(1);
	}
	in.close();
	in.clear();

	in.open(_initial_contig_file_name.c_str());
	if (!in) 
	{
		cerr << _initial_contig_file_name << " : No such file or directory!!!" << endl;
		exit(1);
	}
	in.close();
	in.clear();

}

Gap_Filling::~Gap_Filling()
{
	//
}

//Align two nerghboring contigs
int Gap_Filling::alignment(const string& suffix, const string & prefix)
{

	assert(prefix.size() == suffix.size());

	//int *score = new int[prefix.size()+1];
	//memset(score, 0, (prefix.size()+1)*sizeof(int));
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

	//delete [] score;
	return ret;

}

//compare function
bool myFunc(const string &s1, const string &s2)
{
	return (s1.size() > s2.size());
}

//read contigs and contigs in condensed de Bruijn graph
void Gap_Filling::set_edges()
{
	ifstream infile(edge_file_name.c_str());
	if (!infile)
	{
		cerr << "No such file or directory!!" << endl;
		exit(EXIT_FAILURE);
	}

	string line;
    string seq;
	int cur_index(0);

	//get unique edges for training parameter
#ifdef TRAINING
	vector<string> uniq_edges_for_training;
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
					cout << seq << endl;
					cout << K << endl;
					continue;
				}
#endif

				Edge edge(cur_index, seq.size(), seq, copy_number);

#ifdef TRAINING
				if (copy_number == 1)
				{
					uniq_edges_for_training.push_back(seq);
				}
#endif
				all_cdbg_edges.push_back(edge);
				++cur_index;
				seq.clear();
			}

			fields = Split_String(line).split('\t');
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
	cout << "read contig number = " << all_cdbg_edges.size() << endl;
	infile.close();
	infile.clear();

#ifdef TRAINING
	sort(uniq_edges_for_training.begin(), uniq_edges_for_training.end(),myFunc);
    cout << "Output unique contigs for training..." << endl;

    char buf[10];
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    string fileName;
    fileName += string(buf);
    fileName += string("mer");
    fileName += ".unique_contig_for_training.fasta";

	ofstream out(fileName.c_str());
	for(int i = 0 ; i < MIN_EDGE_COUNT_FOR_TRAINING; ++i)
	{
		out << ">seq_" << i << endl;
		out << uniq_edges_for_training[i] << endl;
	}
	out.close();
	out.clear();
	uniq_edges_for_training.clear();
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

	cout << "all condensed edges in CDBG : " << all_initial_edges.size() << endl;

	infile.close();
	infile.clear();

	cout << "read contigs end..." << endl;
}

//read scaffolds
void Gap_Filling::set_scaffolds()
{

	ifstream infile(scaffold_file_name.c_str());
	if (!infile)
	{
		cerr << scaffold_file_name << " No Such file or directory..." << endl;
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
				temp_split = Split_String(line).split('\t');
				
				if(temp_split.size() == 0)
					continue;

				for(int i = 0 ;i < temp_split.size() ; ++i)
					edges.push_back( atoi(temp_split[i].c_str()));

				edge_number = edges.size();

				getline(infile, line);
				temp_split = Split_String(line).split('\t');

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
	cout << "Input Scaffolds Number = " << all_scaffolds.size() << endl;
}

//run BFS to set current gap 
void Gap_Filling::BFS( const  int gap_index, const  int left_index, const  int right_index, int dis)
{

	vector<vector<int> >  final_sets; //only contain the middle sequence except the two ends
    choice_type prev_sets;
	choice_type next_sets;
	
	dis += (K-1+3*var); //distance constraints
	dis += 40;

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
		if(step > STEP) break;

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

void Gap_Filling::gap_filling()
{
	get_gap_info();
	uniq_candidate_gap_filling();
	multi_candidates_gap_filling();
	get_multi_candidates_gap_local_seq();
}


void Gap_Filling::get_gap_info()
{

	cout << "Begin Gap Filling ... " << endl;
	set_edges();
	set_scaffolds();
	cout << "Set edges and Scaffolds end..." << endl;

	cout << "Begin build condensed de Buijn graph... " << endl;

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
	
    cout << "Begin build the prev and next node of contig in CDBG..." << endl;

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

	cout << "Begin build condensed de Buijn graph... " << endl;

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

    cout << "set edges in CDBG..." << endl;
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
	
    cout << "Build CDBG end..." << endl;

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
	cout << "The number of gaps = " << gap_indexs.size() << endl;
	cout << "Filling the gaps which shares overlap ( < " << overlap << ") = " << pre_gap_info.size() << endl;

	assert(gap_indexs.size() == cur_gap_index);
	assert(gap_indexs.size() == gap_distances.size());

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
		}
		catch(bad_alloc)
		{
			cout << "BFS is too memory-intensive..." << endl;
			fail_gap_info.push_back(i);
		}
	}

	cout << "classify gap end..." << endl;

	cout << "unique candidate gaps number = " << uniq_gap_info.size() << endl;
	cout << "multi candidates gaps number = " << multi_gap_info.size() << endl;
	cout << "failed gaps number = " << fail_gap_info.size() << endl;
	//cout << "gap_state.size() : " << gap_state.size() << endl;
}

//
void Gap_Filling::uniq_candidate_gap_filling()
{
    cout << "Begin filling unique candidates gaps ... " << endl;
	int uniq_gap_number = uniq_gap_info.size();

    string seq;

	for(int i = 0 ;  i < uniq_gap_number ;  ++i )
	{
		int cur_gap_index = uniq_gap_info[i].first;
		vector<int> path = uniq_gap_info[i].second;

		int left_index = gap_indexs[cur_gap_index].first;
		int right_index = gap_indexs[cur_gap_index].second;

        seq.clear();
        //filled by condensed contigs
		if (gap_state.find(cur_gap_index) != gap_state.end())
		{
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

    cout << "Filling unique candidate gaps end..." << endl;
    return;
}

//do not contain the two ends
void Gap_Filling::multi_candidates_gap_filling()
{
    cout << "Begin multi-candidates gap filling... " << endl;

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
			for( int j = 0 ; j < paths.size(); ++j)
			{
				seq.clear();
				seq += all_cdbg_edges[left_index].get_edge_seq();
				
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
			for( int j = 0 ; j < paths.size(); ++j)
			{
				seq.clear();
				seq += all_cdbg_edges[left_index].get_edge_seq();
				
				for( int k = 0; k < paths[j].size(); ++k)
				{
					assert(paths[j][k] < all_initial_edges.size());

					seq += (all_initial_edges[paths[j][k]].get_edge_seq());
					assert(seq.size() >= K);
					seq.resize(seq.size()-(K-1));		
				}
				seq += (all_cdbg_edges[right_index].get_edge_seq());
				seq_s.push_back(seq);
			}
			multi_gap_seq.insert(make_pair(cur_gap_index, seq_s));
		}
	} 
	
	assert(multi_gap_seq.size() == multi_gap_info.size());

	cout << "multi candidate gap filling end ... " << endl;
    
    //clear the initial edges since it would not be used again
	all_initial_edges.clear();
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
		index_path.resize(index_path.size()-1);   //subtract the last '_'
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

/*
 * Filling the unique gap and preOverlap 
 * print the initial gap sequence
 * */

bool myFuncLong(long lhs, long rhs)
{
	return lhs > rhs;
}

bool myFuncPairIntLong(pair<int,long> lhs, pair<int,long> rhs)
{
	return lhs.second > rhs.second;
}	
bool myFuncCompareStr(const string str1, const string str2)
{
	return str1.size() > str2.size();
}

void Gap_Filling::output_initial_scaffolds_seq()
{
	cout << "Output scaffolds sequences..." << endl;

	int index = 0;

	vector<pair<int , long> > initial_scaffolds_len;
	vector<string> initial_scaffolds_seq;

	for(int  i = 0 ; i < all_scaffolds_segments.size() ; ++i)
	{
		for(int j = 0 ; j < all_scaffolds_segments[i].size(); ++j)
		{
			assert(all_scaffolds_segments[i][j].size() > 0);
		}
	}
    string temp_seq;
	
	set<int>::iterator it;
	for(int i = 0 ; i < all_scaffolds_segments.size() ; ++i)
	{
		if(all_scaffolds_segments[i].size() == 1)
		{
			initial_scaffolds_seq.push_back(all_scaffolds_segments[i][0]);
			initial_scaffolds_len.push_back(pair<int,long>(i,all_scaffolds_segments[i][0].size()));
			continue;
		}

		int j = 0;
		long temp_len = 0;

        temp_seq.clear();
		for( ; j < all_scaffolds_segments[i].size()-1; ++j)
		{
		
			temp_seq += all_scaffolds_segments[i][j];
			temp_len += all_scaffolds_segments[i][j].size();
			if (gap_distances[multi_gap[index]] < 0)
			{
				temp_len += 0;
			}
			else
			{
				temp_seq += string(gap_distances[multi_gap[index]], 'N');
				temp_len += gap_distances[multi_gap[index]];
			}

			++index;
		}

		temp_seq += all_scaffolds_segments[i][j];
		temp_len += all_scaffolds_segments[i][j].size();

		initial_scaffolds_len.push_back(pair<int, long>(i,temp_len));
		initial_scaffolds_seq.push_back(temp_seq);

	}
	
	assert(index == multi_gap.size());
	assert(initial_scaffolds_len.size() == all_scaffolds_segments.size());

	sort(initial_scaffolds_len.begin(), initial_scaffolds_len.end(), myFuncPairIntLong);

	char buf[10];
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    string fileName;
    fileName += string(buf);
    fileName += string("mer");
    fileName += ".scaf_len";

    buf[0] = '\0';
    sprintf(buf, "%d",K);
    string fileName2;
    fileName2 += (string(buf));
    fileName2 += "mer";
    fileName2 += ".scaf_seq";
    
	//ofstream outfile((TEMP_FILE_DIR + fileName).c_str());
	//ofstream outfile2((TEMP_FILE_DIR + fileName2).c_str());
	ofstream outfile(fileName.c_str());
	ofstream outfile2(fileName2.c_str());
    outfile << "index\tlength\n" << endl;
    
	for(int i = 0 ; i < initial_scaffolds_len.size(); ++i)
	{
		outfile << initial_scaffolds_len[i].first << "  " << initial_scaffolds_len[i].second << endl;

		outfile2 << ">seq_" << initial_scaffolds_len[i].first << "_" << initial_scaffolds_len[i].second << endl;
		outfile2 << initial_scaffolds_seq[initial_scaffolds_len[i].first] << endl;
	}
	outfile2.close();
	outfile2.clear();
	outfile.close();
	outfile.clear();

	//outfile.open((TEMP_FILE_DIR + "final_result").c_str());
	
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    fileName.clear();
    fileName += string(buf);
    fileName += string("mer");
    fileName += ".final_result";

	outfile.open(fileName.c_str());

	long GENOME = 0;
	for(int i = 0 ; i < initial_scaffolds_len.size(); ++i)
	{
		GENOME += initial_scaffolds_len[i].second;
	}
	
	outfile << "Genome estimated size\t:\t" << GENOME << endl;

	long N50 = GENOME * 0.5 ;
	long N90 = GENOME * 0.9 ;

	outfile << "---Scaffold---"<< endl;

	outfile << "Scaff length : Top 20..." << endl;

	int size = initial_scaffolds_len.size();
	for(int i = 0 ; i < min(20, size); ++i)
	{
		outfile << initial_scaffolds_len[i].second << endl;
	}

	vector<long> sum (initial_scaffolds_len.size(),0);
	
	sum[0] = initial_scaffolds_len[0].second;
	for(int i = 1; i < initial_scaffolds_len.size(); ++i)
	{
		sum[i] +=  sum[i-1] + initial_scaffolds_len[i].second;
	}
	long N50_val = 0;
	long N90_val = 0;

	for(int i = 0 ; i < sum.size(); ++i)
	{
		if(sum[i] >= N50 && N50_val == 0)
		{
			N50_val = initial_scaffolds_len[i].second;
		}
		if (sum[i] >= N90 && N90_val == 0)
		{
			N90_val = initial_scaffolds_len[i].second;
			break;
		}
	}

	outfile << "N50\t:\t" << N50_val << endl;
	outfile << "N90\t:\t" << N90_val << endl;

	outfile << "All gaps number in scaffolds\t:\t" << gap_indexs.size() << endl;
	outfile << "filled gap numbers\t:\t" << uniq_gap_info.size() + multi_gap_info.size() + pre_gap_info.size() << endl;

	outfile << "failed gap numbers\t:\t" << fail_gap_info.size() << endl;
	

	outfile.close();
	outfile.clear();
}

//get multi gap local sequence and merge the short gaps
//but discard the candidate number which is grt than MULTI_CANDIDATE_THRESHOLD
void Gap_Filling::get_multi_candidates_gap_local_seq()
{
	cout << "Get multi candidates gap local seq to PerM and scoring..." << endl;

	get_gap_candidate_number_hash();
    //get all  scaffolds segments and extract the multi gap index to multi_gap
	get_all_scaffolds_segments();
    
	int index = 0;
    
	//store current gap's all local sequences
	vector<string> candidate_seqs;

	//The gaps need to be merged
	vector<int> merged_gap_choices;
	vector<int> merged_gap_indexs;

    set<int> merge_segment_index; //record the merged segment index

	assert(new_gap_indexs.size() == 0);

    char buf[10];
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    string fileName;
    fileName += string(buf);
    fileName += "mer";
    fileName += ".multi_candidates_seq_to_PerM.fasta";


	//outfile.open((TEMP_FILE_DIR + fileName).c_str());
	ofstream out_multi_gaps(fileName.c_str());
	if (!out_multi_gaps)
	{
		throw runtime_error("No such file or directory!!!");
		exit(1);
	}
	
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    fileName.clear();
    fileName += string(buf);
    fileName += "mer";
    fileName += ".multi_gap_count_to_score";

	//ofstream outfile2((TEMP_FILE_DIR + fileName).c_str());
	ofstream out_gap_count(fileName.c_str());
	if (!outfile2)
	{
		throw runtime_error("No such file or directory!!!");
		exit(1);
	}


    buf[0] = '\0';
    sprintf(buf, "%d",K);
    fileName.clear();
    fileName += string(buf);
    fileName += "mer";
    fileName += ".multi_scaff_segment_to_score";
	ofstream out_scaff_seg(fileName.c_str());
    
    srand(time(NULL));

	int threshold = MU + var*3;
    int cur_gap_index = 0 ; //record the gap index of current gap pair

	vector<string> cur_scaffold_seq;

	for(int i = 0; i < all_scaffolds_segments.size(); ++i)
	{
		vector<string> seqs = all_scaffolds_segments[i];
		assert(seqs.size() > 0);

		cur_scaffold_seq.clear();

		if (seqs.size() == 1)
		{
			cur_scaffold_seq.push_back(seqs[0]);
		//	new_all_scaffolds_segments.push_back(cur_scaffold_seq);
            out_scaff_seg << ">seq " << i << endl;
            out_scaff_seg << seqs[0] << endl;
			continue;
		}
		int j = 1, cur = 0 ;

		cur_scaffold_seq.push_back(seqs[cur]);

		int merge_gap_count = 0;
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

			vector<string> adjacent;

			for(int k = cur; k < j; ++k)
			{
				merged_gap_indexs.push_back(multi_gap[index]);
				merged_gap_choices.push_back((multi_gap_seq[multi_gap[index]]).size());
				adjacent.push_back(seqs[k]);
                ++index;
			}

			adjacent.push_back(seqs[j]);
			assert(adjacent.size() == (j-cur+1));
			assert(adjacent.size() == merged_gap_indexs.size()+1);
			assert(merged_gap_indexs.size() == (j-cur));
			assert(merged_gap_choices.size() == (j-cur));
			
			//Record the index sequence,as 0_0_1
			vector<string> indexs;
			string index_path;
			string path;

            //subtract the suffix sequence of MU+3*var
			if (adjacent[0].size() > threshold)
			{
				adjacent[0] = adjacent[0].substr(adjacent[0].size()-threshold);
			}
            //subtract the prefix sequence of MU+3*var
			if (adjacent[adjacent.size()-1].size() > threshold)
			{
				adjacent[adjacent.size()-1].resize(threshold);
			}
	
			DFS(merged_gap_indexs, merged_gap_choices, indexs, 0, candidate_seqs, index_path, path, adjacent);

			new_gap_indexs.push_back(make_pair(
						gap_indexs[merged_gap_indexs[0]].first , 
						gap_indexs[merged_gap_indexs[merged_gap_indexs.size()-1]].second));

			for(int l = 0 ; l < merged_gap_indexs.size(); ++l)
			{
				assert(multi_gap_candidate_seq.find(merged_gap_indexs[l]) == multi_gap_candidate_seq.end());
				multi_gap_candidate_seq[merged_gap_indexs[l]] = candidate_seqs;
			}
			new_multi_gap_candidate_seq.insert(make_pair(new_gap_indexs.size()-1, candidate_seqs));

			candidate_seqs.clear();
			cur = j;
			j++;

			cur_scaffold_seq.push_back(seqs[cur]);
		}
		//new_all_scaffolds_segments.push_back(cur_scaffold_seq);
        outfile << ">seq  " << i << endl;
        for(int j=0; j < cur_scaffold_seq.size(); ++j)
        {
          if(merge.find(j) != merge.end())
          {
            out_scaff_seg << cur_scaffold_seq[j];
            out_scaff_seg << new_multi_gap_candidates_seq[index]
          }
        }
	}
		
//	cout << "new_all_scaffolds_segments end..." << endl;
//	cout << "new scaffolds size : " << new_all_scaffolds_segments.size() << endl;

//	assert(new_all_scaffolds_segments.size() == all_scaffolds_segments.size());
//	assert(new_all_scaffolds_segments.size() == all_scaffolds.size());
	assert(new_multi_gap_candidate_seq.size() == new_gap_indexs.size());
	assert(index == multi_gap.size());

	//reverse new gap pair index to hash
	for(int i = 0 ; i < new_gap_indexs.size(); ++i)
	{
		assert(new_gap_indexs_hash.find(new_gap_indexs[i]) == new_gap_indexs_hash.end());
		new_gap_indexs_hash.insert(make_pair(new_gap_indexs[i], i));
	}
	assert(new_gap_indexs_hash.size() == new_gap_indexs.size());
/*
	cout << "output multi candidates gaps sequences to perM and score" << endl;

    char buf[10];
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    string fileName;
    fileName += string(buf);
    fileName += "mer";
    fileName += ".multi_candidates_seq_to_PerM.fasta";


	//outfile.open((TEMP_FILE_DIR + fileName).c_str());
	ofstream outfile(fileName.c_str());
	if (!outfile)
	{
		throw runtime_error("No such file or directory!!!");
		exit(1);
	}
	
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    fileName.clear();
    fileName += string(buf);
    fileName += "mer";
    fileName += ".multi_gap_count_to_score";
	
	//ofstream outfile2((TEMP_FILE_DIR + fileName).c_str());
	ofstream outfile2(fileName.c_str());
	if (!outfile2)
	{
		throw runtime_error("No such file or directory!!!");
		exit(1);
	}

	index = 0;

    int temp_count = 0;
    set<int > merge_segment_index; // current segment[i][j] need to be merged t its previous

	for(int i = 0 ; i < new_multi_gap_candidate_seq.size(); ++i)
	{
		temp_count = new_multi_gap_candidate_seq[i].size();
        if(temp_count > MAX_CHOICE)
        {
          merge_segment_index.insert(i);
          continue;
        }
		outfile2 << temp_count << endl;

		for(int j = 0 ; j < new_multi_gap_candidate_seq[i].size(); ++j)
		{
			outfile << ">seq_" << index << "_";
			outfile << j << endl;
			outfile<< new_multi_gap_candidate_seq[i][j] << endl;
		}
		++index;
	}
	outfile.close();
	outfile.clear();

	outfile2.close();
	outfile2.clear();
    */

	//get initial scaffolds and no multi gaps and failed gaps which is filled by 'N'
	output_initial_scaffolds_seq();
    //cout << "!!! new multi gaps : " << new_multi_gap_candidate_seq.size() << endl;
 //   merge_segments(merge_segment_index);
	//cout << "!!!need to merge : " << merge_segment_index.size() << endl;
}

void Gap_Filling::merge_segments(const set<int > &merge) 
{
    char buf[10];
    buf[0] = '\0';
    sprintf(buf, "%d",K);
    string fileName;
    fileName += string(buf);
    fileName += "mer";
    fileName += ".multi_scaff_segment_to_score";

	//ofstream outfile((TEMP_FILE_DIR+fileName).c_str());
	ofstream outfile(fileName.c_str());
    
    srand(time(NULL));

    int count = 0;
	for(int i = 0 ; i < new_all_scaffolds_segments.size() ; ++i)
	{
        count += (new_all_scaffolds_segments[i].size()-1);
	}
	
    //cout << "!!!!count : " << count << endl;

    int index = 0;
    int index2 = 0;
	for(int i = 0 ; i < new_all_scaffolds_segments.size() ; ++i)
	{
		outfile << ">seq" << i << endl;
		for(int j = 0 ; j < new_all_scaffolds_segments[i].size() ; ++j)
		{
          if(merge.find(index) != merge.end())
          {
              outfile << new_all_scaffolds_segments[i][j];
              outfile << new_multi_gap_candidate_seq[index][rand()%(new_multi_gap_candidate_seq[index].size())];
              ++index;
          }
		  else	
          {
            outfile << new_all_scaffolds_segments[i][j] << endl;
            ++index;
            ++index2;
          }
		}
	}
	outfile.close();
	outfile.clear();
    //cout << "!!!!!index : " << index << endl;
    //cout << "!!!!!index2 : " << index2-new_all_scaffolds_segments.size() << endl;
    //cout << "!!!!!index : " << index-new_all_scaffolds_segments.size() << endl;
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

//scaffolds_seq with several segments since the multi gap 
//get scaffold segments 
void Gap_Filling::get_all_scaffolds_segments()
{

	int cur_gap_index = 0;
	string seq;
	vector<string> scaffold_segment;
	map<pair<int,int>,int>::const_iterator it;

	for(int i = 0; i != all_scaffolds.size(); ++i)
	{

		seq.clear();
		scaffold_segment.clear();

		vector<int> edge_indexs = all_scaffolds[i].get_edges_of_scaffold();
		assert(edge_indexs.size() > 0);

		seq += all_cdbg_edges[edge_indexs[0]].get_edge_seq();

		if (edge_indexs.size() == 1)
		{
			scaffold_segment.push_back(seq);
			assert(scaffold_segment.size() > 0);
			all_scaffolds_segments.push_back(scaffold_segment);
			continue;
		}

		for(int j = 0; j+1 < edge_indexs.size(); ++j)
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
				assert(seq.size() > abs(it->second));
				seq.resize(seq.size() + it->second);
				seq += all_cdbg_edges[right_index].get_edge_seq();
			}
			//unique
			else if(it->second == 1 )
			{
				assert(seq.size() >= left_len);
				seq.resize(seq.size() - left_len);
				seq += uniq_gap_seq[cur_gap_index];
			}
			//failed gap
			else if(it->second == 0)
			{
              //for failed gap , the count of 'N' is determined by the distance estimation
				if(gap_distances[cur_gap_index] > (MU+3*var))
				{
					seq += string(MU+3*var,'N');
				}
				else if (gap_distances[cur_gap_index] < 0)
				{
					;
				}
				else
				{
					seq += string(gap_distances[cur_gap_index],'N');
				}
				seq += all_cdbg_edges[right_index].get_edge_seq();
			}
			//multi_gap
			else
			{
				assert(it->second <= MAX_CHOICE);
				scaffold_segment.push_back(seq);
				multi_gap.push_back(cur_gap_index);
				seq.clear();
				seq += all_cdbg_edges[right_index].get_edge_seq();
			}
			cur_gap_index++;
		}
		scaffold_segment.push_back(seq);
		assert(scaffold_segment.size() > 0);
		all_scaffolds_segments.push_back(scaffold_segment);
	}
	assert(all_scaffolds_segments.size() == all_scaffolds.size());
	assert(cur_gap_index == gap_indexs.size());
	assert(multi_gap.size() == multi_gap_info.size());
    //cout << "!!!!! multi_gap.size() : ( = 611)" << multi_gap.size() << endl;
}


//operator<<
ostream &operator<< (ostream &os, const Gap_Filling &obj )
{
	os << "pre_gap_info..." << endl;
	for( map< int , int>::const_iterator it = (obj.pre_gap_info).begin();
			it != (obj.pre_gap_info).end();
			++it)
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
		copy( (it->second).begin(),
				(it->second).end(),
				ostream_iterator<string>(os, "\n"));
		os << endl;
	}
	os << "each scaffold segment..." << endl;

	int sum = 0 ;
	int index = 0;
	assert(obj.all_scaffolds_segments.size() == obj.all_scaffolds.size());
	for(int i = 0; i < obj.all_scaffolds_segments.size() ; ++i)
	{
		assert(obj.all_scaffolds_segments[i].size() > 0);

		//cout << obj.all_scaffolds_segments[i].size () << endl;
		if (obj.all_scaffolds_segments[i].size() > 1)
		{
			sum += (obj.all_scaffolds_segments[i].size()-1);
		}
		os << "i = " << i << endl;

		int j = 0;
		for(; j < obj.all_scaffolds_segments[i].size()-1; ++j)
		{
			//cout << "current pair = " << obj.gap_indexs[obj.multi_gap[j+index]].first << " " << obj.gap_indexs[obj.multi_gap[j+index]].second << endl;
		}
		index += j;
		/*
		//cout << "DEBUG3 , i = " << i << endl;
		for(int l = 0 ; l < obj.all_scaffolds_segments[i].size(); ++l) 
		{
			cout << obj.all_scaffolds_segments[i][l].size() << " ";
		}
		cout << endl;
		*/
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
			//			cout  
			//				<< "current gap = " 
			//				<< (obj.gap_indexs)[i].first
			//				<< " " 
			//				<< (obj.gap_indexs)[i].second 
			//				<< endl;	
			os << obj.multi_gap[index] << "\t";
		}
	}

	cout << "print each scaffold segments end..." << endl;
	return os;

}
