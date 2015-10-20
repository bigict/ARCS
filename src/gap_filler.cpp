#include "gap_filler.h"
#include "component.h"
#include "condensed_debruijn_graph.h"
#include "condensed_debruijn_graph_reader.h"
#include "constant.h"
#include "contigs.h"
#include "kmer.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.GapFiller"));

#define TRAINING
#define RESULT

//gap filling
void GapFiller::fill() {
	get_gap_info();
	uniq_candidate_gap_filling();
	multi_candidates_gap_filling();
	output_initial_scaffolds_seq();
	get_multi_candidates_gap_local_seq();
}

//Align two nerghboring contigs
size_t GapFiller::alignment(const string& suffix, const string& prefix) {
	assert(prefix.size() == suffix.size());

    std::vector< size_t > score(prefix.size() + 1, 0);

	for (size_t i = 1; i <= suffix.size(); ++i) {
		for (size_t j = prefix.size(); j > 0; --j) {
			if (prefix[j - 1] == suffix[i - 1]) {
				score[j] = score[j - 1] + 1;
			} else {
				score[j] = 0;
			}
		}
	}

	size_t ret = 0;
	for (size_t j = 0; j <= prefix.size(); ++j) {
		if (score[j] == j) {
			if (ret < score[j]) {
				ret = score[j];
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

bool GapFiller::input_scaffold(const std::string& file) {
    std::ifstream stream(file.c_str());
    if (!stream) {
        LOG4CXX_WARN(logger, boost::format("%s: No such file or directory!!!") % file);
        return false;
    }

    ComponentReader reader(stream);
    Component component;
    while (reader.read(component)) {
        if (!component.gaps.empty()) {
            component.gaps.pop_back();
        }
        //Scaffold scaff(curIdx, (int)component.contigs.size(), component.contigs, component.gaps);
        all_scaffolds.push_back(component);
    }
    return true;
}

bool GapFiller::input_contigs(const std::string& file) {
	//get unique edges for training parameter
#ifdef TRAINING
    priority_queue<string, vector<string>, mycomparison> unique_contigs_for_training;
#endif

    // cdbg_copy_number.fa
    {
        std::ifstream stream(file.c_str());
        if (!stream) {
            LOG4CXX_WARN(logger, boost::format("%s: No such file or directory!!!") % file);
            return false;
        }
        size_t curIdx = 0;
        ContigReader reader(stream);
        Contig contig;
        while (reader.read(contig)) {
#ifdef TRAINING
            if (contig.copy_num == 1) {
                unique_contigs_for_training.push(contig.seq);
            }
#endif
            _uniq_graph.addEdge(contig.seq, contig.copy_num);
            ++curIdx;
        }

	    LOG4CXX_DEBUG(logger, boost::format("all contig number=%d") % _uniq_graph._indexer.size());
    }

#ifdef TRAINING
    {
	    LOG4CXX_DEBUG(logger, "Output unique contigs for training...");
        std::string file = boost::str(boost::format("%dmer.unique_contig_for_training.fasta") % K);
        std::ofstream stream(file.c_str());
        if (!stream) {
            LOG4CXX_WARN(logger, boost::format("%s: No such file or directory!!!") % file);
            return false;
        }
        for (size_t i = 0; i < MIN_EDGE_COUNT_FOR_TRAINING && !unique_contigs_for_training.empty(); ++i) {
            stream << boost::format(">seq_%d") % i << std::endl;
            stream << unique_contigs_for_training.top() << std::endl;
            unique_contigs_for_training.pop();
        }
        while (!unique_contigs_for_training.empty()) {
            unique_contigs_for_training.pop();
        }
    }
#endif
    return true;
}

bool GapFiller::input_debruijn(const std::string& file) {
    // condensed_de_bruijn_graph_before_trimming.data
    std::ifstream stream(file.c_str());
    if (!stream) {
        LOG4CXX_WARN(logger, boost::format("%s: No such file or directory!!!") % file);
        return false;
    }
    CondensedDeBruijnGraphReader reader(stream);
    CondensedDeBruijnGraphEdge debruijn;
    while (reader.read(debruijn)) {
        _all_graph.addEdge(debruijn.seq, 0);
    }
    LOG4CXX_DEBUG(logger, boost::format("All condensed edges number=%d") % _all_graph._indexer.size());
}

void GapFiller::BFS(const size_t i, const size_t j, std::vector< std::vector< size_t > >& pathlist) {
}
// Run BFS to fill current gap 
void GapFiller::BFS(const int gap_index, const  int left_index, const  int right_index, int dis) {

	vector<vector<int> >  final_sets; //only contain the middle sequence except two ends
	typedef std::vector< std::pair< std::vector< int >, int > > ChoiceList;
	ChoiceList prev_sets;
	ChoiceList next_sets;

	dis += (K-1+3*var); //distance constraints
	dis += 30;

	int step = 0;

    {
        const CondensedDeBruijnGraph::CondensedEdge& edge = _uniq_graph._indexer[left_index];
        std::string suffix = edge.seq.substr(edge.seq.length() - (K - 1), K - 1);
        CondensedDeBruijnGraph::NodeList::const_iterator i = _uniq_graph._parents.find(suffix);
        if (i != _uniq_graph._parents.end()) {
            for (CondensedDeBruijnGraph::EdgeList::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
                size_t length = _uniq_graph._indexer[j->second].seq.length() - (K - 1);
                if (j->second == right_index) {
                    final_sets.push_back(vector<int>(1,j->second));
                }
                prev_sets.push_back(make_pair( vector< int > ( 1 , j->second) , length )) ;
            }
        }
        ++step;
    }

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

		for (int k = 0 ; k < prev_sets.size(); ++k)
		{
			vector<int> path = prev_sets[k].first;
			int index = path[path.size()-1];
			int len = prev_sets[k].second;
			if ( len > dis ) continue;

            const CondensedDeBruijnGraph::CondensedEdge& edge = _uniq_graph._indexer[index];
            std::string suffix = edge.seq.substr(edge.seq.length() - (K - 1), K - 1);
            char prev = edge.seq[edge.seq.length() - K];
            CondensedDeBruijnGraph::NodeList::const_iterator i = _uniq_graph._parents.find(suffix);
            if (i != _uniq_graph._parents.end()) {
                for (CondensedDeBruijnGraph::EdgeList::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
                    if (j->second == right_index) {
                        final_sets.push_back(path);
                    }

                    int len2 = len + _uniq_graph._indexer[j->second].seq.length() - (K - 1);
                    path.push_back(j->second);
                    next_sets.push_back(make_pair(path,len2));
                    path.pop_back();
                }
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

    const std::string& lsequence = _uniq_graph._indexer[left_index].seq;
    const std::string& rsequence = _uniq_graph._indexer[right_index].seq;
    std::string start = lsequence.substr(lsequence.length() - (K - 1), K - 1);
    std::string end = rsequence.substr(0, K - 1);

	set<int> ends;

    {
        CondensedDeBruijnGraph::NodeList::const_iterator i = _all_graph._children.find(end);
        if (i == _all_graph._children.end()) {
		    fail_gap_info.push_back(gap_index);
            return;
        }
        for (CondensedDeBruijnGraph::EdgeList::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
            ends.insert(j->second);
        }
    }
	if (ends.empty()) {
		fail_gap_info.push_back(gap_index);
		return;
	}

	//find the first step indexs in initial
    {
        CondensedDeBruijnGraph::NodeList::const_iterator i = _all_graph._parents.find(start);
        if (i == _all_graph._parents.end()) {
            fail_gap_info.push_back(gap_index);
            prev_sets.clear();
            next_sets.clear();
            final_sets.clear();

            return;
        }
        for (CondensedDeBruijnGraph::EdgeList::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
            int len = _all_graph._indexer[j->second].seq.length() - (K-1);

            if ( ends.find(j->second) != ends.end() )
            {
                final_sets.push_back(vector<int>(1, j->second));
            }
            prev_sets.push_back(make_pair( vector< int > ( 1 , j->second ) , len )) ;
        }
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
		for (int k = 0; k < prev_sets.size(); ++k)
		{
			vector<int> path = prev_sets[k].first;
			int index = path[path.size()-1];
			int len = prev_sets[k].second;

			if(len > dis)
			{
				continue;
			}

            const CondensedDeBruijnGraph::CondensedEdge& edge = _all_graph._indexer[index];
            std::string suffix = edge.seq.substr(edge.seq.length() - (K - 1), K - 1);
            char prev = edge.seq[edge.seq.length() - K];
            CondensedDeBruijnGraph::NodeList::const_iterator i = _all_graph._parents.find(suffix);
            if (i != _all_graph._parents.end()) {
                for (CondensedDeBruijnGraph::EdgeList::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
                    if (ends.find(j->second) != ends.end()) {
                        final_sets.push_back(path);
                    }
                    int len2 = len + _all_graph._indexer[j->second].seq.length()-(K-1);
                    path.push_back(j->second);
                    next_sets.push_back(make_pair(path,len2));
                    path.pop_back();
                }
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

//get gap info and classify gaps
void GapFiller::get_gap_info() {
    LOG4CXX_DEBUG(logger, "Begin Gap Filling ... ");

    std::set< std::pair< int, int > > visited; //record the processed gaps
	int cur_gap_index(0);

	for (size_t i = 0; i < all_scaffolds.size(); ++i) {
		if (all_scaffolds[i].edge_num() == 1) { // no gaps
			continue;
		}
        Component::ContigIdList edge_indexs = all_scaffolds[i].contigs;

		for (size_t j = 1; j < edge_indexs.size(); ++j) {
			int left_index = edge_indexs[j - 1];
			int right_index = edge_indexs[j];

			gap_indexs.push_back(make_pair(left_index, right_index));
			gap_distances.push_back(all_scaffolds[i].gaps[j - 1]);

            const std::string& lsequence = _uniq_graph._indexer[left_index].seq;
            const std::string& rsequence = _uniq_graph._indexer[right_index].seq;
            std::string suffix = lsequence.substr(lsequence.length() - (K - 1), K - 1);
            std::string prefix = rsequence.substr(0, K - 1);

			int temp_overlap = 0;
			if ((temp_overlap = alignment(suffix, prefix)) >= overlap)
			{
				assert(visited.find(make_pair(left_index, right_index)) == visited.end());
				visited.insert(make_pair(left_index, right_index));
				pre_gap_info.insert(make_pair(cur_gap_index,temp_overlap));
			}
			cur_gap_index++;
		}
	}
	LOG4CXX_DEBUG(logger, boost::format("The number of gaps = %d") % gap_indexs.size());
	LOG4CXX_DEBUG(logger, boost::format("Filling the gaps which shares overlap ( < %d) = %d") % overlap % pre_gap_info.size());

	assert(gap_indexs.size() == cur_gap_index);
	assert(gap_indexs.size() == gap_distances.size());

#ifdef RESULT 
    int count_base_in_gap=0;
#endif

	// for each gap , run BFS and get uniq_gap_info and multi_gap_info
	for (size_t i = 0 ; i < gap_indexs.size(); ++i) {

		int left_index = gap_indexs[i].first;
		int right_index = gap_indexs[i].second;

		if (visited.find(gap_indexs[i]) != visited.end()) {
			continue;
		}

		int distance = gap_distances[i];
		if (distance > (MU + 3*var )) {
			fail_gap_info.push_back(i);
			continue;
		}

		try {
			BFS(i, left_index, right_index, distance);
#ifdef RESULT 
            count_base_in_gap += distance;
#endif
		} catch(bad_alloc) {
			cerr << "[Info] BFS is too memory-intensive, ignoring..." << endl;
			fail_gap_info.push_back(i);
		}
	}

#ifdef  RESULT
    LOG4CXX_DEBUG(logger, boost::format("Base counts in gaps : %d") % count_base_in_gap);
#endif

	LOG4CXX_DEBUG(logger, boost::format("Unique candidate gaps number = %d") % uniq_gap_info.size());
	LOG4CXX_DEBUG(logger, boost::format("Multi candidates gaps number = %d") % multi_gap_info.size());
	LOG4CXX_DEBUG(logger, boost::format("Failed gaps number = %d") % fail_gap_info.size());
}

//Filling the unique candidate gaps
//contain two ends
void GapFiller::uniq_candidate_gap_filling() {
	LOG4CXX_DEBUG(logger, "Start filling unique candidates gaps ... ");

	int uniq_gap_number = uniq_gap_info.size();

	string seq;

	for (size_t i = 0;  i < uniq_gap_info.size(); ++i) {
		size_t gap_index = uniq_gap_info[i].first;
        std::vector< int >& path = uniq_gap_info[i].second;

		int left_index = gap_indexs[gap_index].first;
		int right_index = gap_indexs[gap_index].second;

        std::string seq = _uniq_graph._indexer[left_index].seq;
        assert(seq.size() >= K);
        seq.resize(seq.length() - (K - 1));

        CondensedDeBruijnGraph* graph = NULL;
		if (gap_state.find(gap_index) != gap_state.end()) { //filled by condensed contigs
            graph = &_uniq_graph;
		} else { //filled by initial contigs
            graph = &_all_graph;
		}
        for (size_t j = 0; j < path.size(); ++j) {
            assert(graph->_indexer[path[j]].seq.length() >= K);
            seq += graph->_indexer[path[j]].seq;
            assert(seq.size() >= K);
            seq.resize(seq.size() - (K - 1));
        }

        seq += _uniq_graph._indexer[right_index].seq;
        uniq_gap_seq.insert(std::make_pair(gap_index, seq));
	}
	assert(uniq_gap_seq.size() == uniq_gap_info.size());
	LOG4CXX_DEBUG(logger, "Fill unique candidate gaps end...");
}

//Do not contain the two ends
void GapFiller::multi_candidates_gap_filling() {
	LOG4CXX_DEBUG(logger, "Start multi candidates gap filling... ");

	for(size_t i = 0; i < multi_gap_info.size(); ++i) {
		assert(multi_gap_info[i].second.size() <= MAX_CHOICE);

		int gap_index = multi_gap_info[i].first;
		int left_index = gap_indexs[gap_index].first;
		int right_index = gap_indexs[gap_index].second;

        std::vector<vector<int> > & paths = multi_gap_info[i].second;
		assert(paths.size() > 1);

        CondensedDeBruijnGraph* graph = NULL;

		if (gap_state.find(gap_index) != gap_state.end()) { //The gap is filled by contigs in cdbg 
            graph = &_uniq_graph;
		} else { //The gap is filled by contigs in initial CDBG
            graph = &_all_graph;
        }

        std::vector< std::string > seqlist;

        for (size_t j = 0 ; j < paths.size(); ++j) {
            std::string seq = graph->_indexer[left_index].seq;
            assert(seq.size() >= K);
            seq.resize(seq.size() - (K - 1));

            for (size_t k = 0; k < paths[j].size(); ++k) {
                seq += graph->_indexer[paths[j][k]].seq;
                assert(seq.size() >= K);
                seq.resize(seq.size() - (K - 1));
            }
            seq += graph->_indexer[right_index].seq;

            seqlist.push_back(seq);
        }

        multi_gap_seq.insert(make_pair(gap_index, seqlist));
	} 

	assert(multi_gap_seq.size() == multi_gap_info.size());

	LOG4CXX_DEBUG(logger, "Multi candidate gap filling end ... ");
}

/*
 * build hash table  map<pair<int, int> , int >
 * uniq -- 1
 * multi --  >1
 * pre -- -(overlap)
 * fail -- 0
 * */

void GapFiller::get_gap_candidate_number_hash() {

	//unique gap -- 1
	for (size_t i = 0; i < uniq_gap_info.size(); ++i) {
        std::pair<int, int> edge_pair = gap_indexs[uniq_gap_info[i].first];
		gap_candidate_number_hash.insert(std::make_pair(edge_pair, 1));
	}
	assert(gap_candidate_number_hash.size() == uniq_gap_info.size());

	//multi gap -- #choices
	for (size_t i = 0; i < multi_gap_info.size(); ++i) {
        std::pair<int, int> edge_pair = gap_indexs[multi_gap_info[i].first];
		int candidate_number = multi_gap_info[i].second.size();

		assert(candidate_number <= MAX_CHOICE);
		gap_candidate_number_hash.insert(std::make_pair(edge_pair, candidate_number));
	}
	assert(gap_candidate_number_hash.size() == uniq_gap_info.size()+ multi_gap_info.size());

	//pre gap info -- -(overlap)
	for (std::map<int,int>::iterator it = pre_gap_info.begin(); it != pre_gap_info.end(); ++it) {
        std::pair<int,int> edge_pair = gap_indexs[it->first];
		gap_candidate_number_hash.insert(std::make_pair(edge_pair, -(it->second)));
	}
	//fail gap info -- 0
	for (size_t i = 0; i < fail_gap_info.size(); ++i) {
        std::pair<int, int> edge_pair = gap_indexs[fail_gap_info[i]];
		gap_candidate_number_hash.insert(std::make_pair(edge_pair, 0));
	}

	assert(gap_candidate_number_hash.size() == uniq_gap_info.size()+ multi_gap_info.size() + pre_gap_info.size() + fail_gap_info.size());
	assert(gap_candidate_number_hash.size() == gap_indexs.size());
}

bool myFuncPairIntLong(pair<int,long> lhs, pair<int,long> rhs)
{
	return lhs.second > rhs.second;
}

//get initial scaffolds with no filling gap filling and using 'N' to fill failed gaps
void GapFiller::output_initial_scaffolds_seq() {
	//get gap number hash
	get_gap_candidate_number_hash();

	int cur_gap_index = 0 ;

	long GENOME = 0;		//estimate the genome size using total length of scaffold
    std::string file = boost::str(boost::format("%dmer.scaf_seq_with_gaps") % K);

    //(index , length) pair of scaffold
    vector<pair<int, long> > scaf_index_len ;
    ofstream out_scaf_seq(file.c_str());
    if (!out_scaf_seq) {
        cerr << "[Info] Create " << file << " failed..." << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < all_scaffolds.size(); ++i) {
        std::string seq;         //scaffold sequence
        std::string one_segment; //one segment of scaffold since the multi gaps

        //for scaffolds segments
        std::vector< std::string > cur_scaff_segments;

        std::vector< size_t > edge_indexs = all_scaffolds[i].contigs;
        assert(edge_indexs.size() > 0);
        seq += _uniq_graph._indexer[edge_indexs[0]].seq;
        one_segment += _uniq_graph._indexer[edge_indexs[0]].seq;

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

            int left_len = _uniq_graph._indexer[left_index].seq.length();
            int right_len = _uniq_graph._indexer[right_index].seq.length();

            pair<int,int> edge_pair = make_pair(left_index,right_index);

            assert(gap_indexs[cur_gap_index].first == left_index && gap_indexs[cur_gap_index].second == right_index);

            map<pair<int,int>, int>::iterator it = gap_candidate_number_hash.find(edge_pair);

            assert( it != gap_candidate_number_hash.end());
            //overlap
            if(it->second < 0)
            {
                assert(one_segment.size() > abs(it->second));
                one_segment.resize(one_segment.size() + it->second);
                one_segment += _uniq_graph._indexer[right_index].seq;
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
                one_segment += _uniq_graph._indexer[right_index].seq;
            }
            //multi_gap
            else
            {
                assert(it->second <= MAX_CHOICE);
                cur_scaff_segments.push_back(one_segment);
                initial_multi_gaps.push_back(cur_gap_index);
                seq += one_segment;
                one_segment.clear();
                one_segment += _uniq_graph._indexer[right_index].seq;
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
    std::string fileName = boost::str(boost::format("%dmer.scaf_len") % K);

	ofstream out_scaf_len(fileName.c_str());
	if(!out_scaf_len)
	{
		cerr << "[Info] " << fileName << " No such file or directory!!!" << endl;
		exit(EXIT_FAILURE);
	}

    fileName = boost::str(boost::format("%dmer.final_result") % K);

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
}

//DFS to get the indexs of new multi gaps to vetify the corresponding index of initial gaps
void GapFiller::DFS(
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
void GapFiller::get_multi_candidates_gap_local_seq() {
    //DEBUG
    {
        int debug_count = 0;
        for(int i=0;i<all_scaffolds_segments.size() ;++i)
        {
            debug_count += (all_scaffolds_segments[i].size()-1);
        }
        assert(debug_count == multi_gap_info.size());
    }

	//clear useless data structure
	uniq_gap_info.clear();
	uniq_gap_seq.clear();

    //cout << all_scaffolds_segments.size() << endl;

	LOG4CXX_DEBUG(logger, "Output multi candidates gap local seq to scoring...");

	string fileName;

    fileName = boost::str(boost::format("%dmer.multi_candidates_seq_to_PerM.fasta") % K);
	ofstream out_multi_gaps(fileName.c_str());
	if (!out_multi_gaps)
	{
		cerr << "[Info] Create " << fileName << " No such file or directory!!!" << endl;
		exit(EXIT_FAILURE);
	}

    fileName = boost::str(boost::format("%dmer.multi_gap_count_to_score") % K);
	ofstream out_gap_count(fileName.c_str());
	if (!out_gap_count)
	{
		cerr << "[Info] " << fileName << " No such file or directory!!!" << endl;
		exit(EXIT_FAILURE);
	}

	fileName = boost::str(boost::format("%dmer.multi_scaff_segment_to_score") % K);
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
std::ostream& operator<<(std::ostream& os, const GapFiller &obj) {
	os << "Pre_gap_info..." << std::endl;
	for (std::map< int, int >::const_iterator it = obj.pre_gap_info.begin(); it != obj.pre_gap_info.end(); ++it) {
        os << boost::format("(%d,%d)%d") % obj.gap_indexs[it->first].first % obj.gap_indexs[it->first].second % it->second << std::endl;
	}

	os << "uniq_gap_info..." << std::endl;
	for (size_t  i = 0; i < obj.uniq_gap_info.size(); ++i) {
        os << boost::format("%d %d : ") % obj.gap_indexs[obj.uniq_gap_info[i].first].first % obj.gap_indexs[obj.uniq_gap_info[i].first].second;
        std::copy(
                obj.uniq_gap_info[i].second.begin(), 
                obj.uniq_gap_info[i].second.end(), 
                std::ostream_iterator< int >(os, "\t"));
		os << std::endl;
	}

	os << "multi_gap_info..." << std::endl;
	for (size_t i = 0 ; i < obj.multi_gap_info.size(); ++i) {
        os << boost::format("( %d,%d )  : ") % obj.gap_indexs[obj.multi_gap_info[i].first].first % obj.gap_indexs[obj.multi_gap_info[i].first].second;
		for (size_t j = 0 ; j < (obj.multi_gap_info[i].second).size(); ++j) {
            os << boost::format("%d = ") % j;
            std::copy(
                    obj.multi_gap_info[i].second[j].begin(), 
                    obj.multi_gap_info[i].second[j].end(), 
                    std::ostream_iterator< int > (os, "\t")
                    );
			os << std::endl;
		}
	}

	os << "uniq_gap_seq\n " << std::endl;
	for (std::map< int, std::string >::const_iterator it = obj.uniq_gap_seq.begin();
			it != obj.uniq_gap_seq.end();
			++it) {
        os << boost::format("( %d,%d )") % obj.gap_indexs[it->first].first % obj.gap_indexs[it->first].second << std::endl;
        os << it->second << std::endl;
	}

	os << "multi_gap_seq\n" << std::endl;
	for (std::map< int, std::vector< std::string > >::const_iterator it = obj.multi_gap_seq.begin();
			it != obj.multi_gap_seq.end();
			++it) {

        os << boost::format("( %d,%d )") % obj.gap_indexs[it->first].first % obj.gap_indexs[it->first].second << std::endl;
        std::copy(
                it->second.begin(), 
                it->second.end(), 
                std::ostream_iterator< std::string >(os, "\n")
                );
		os << std::endl;
	}

	os << "each scaffold segment..." << std::endl;
	assert(obj.all_scaffolds_segments.size() == obj.all_scaffolds.size());
	for (size_t i = 0; i < obj.all_scaffolds_segments.size() ; ++i) {
		assert(obj.all_scaffolds_segments[i].size() > 0);

        os << boost::format("i = %d") % i << std::endl;

        std::copy(
                obj.all_scaffolds_segments[i].begin(), 
                obj.all_scaffolds_segments[i].end(),
                std::ostream_iterator< std::string >(os, "\n")
                );
		os << endl;
	}

	for (std::vector< std::vector< std::string > >::const_iterator 
			it = obj.all_scaffolds_segments.begin();
			it != obj.all_scaffolds_segments.end();
			++it) {

        std::copy(it->begin(), it->end(),  std::ostream_iterator< std::string >(os, "\n"));
		for (size_t i = 0; i < it->size()-1; ++i) {
			os 
				<< "current gap = " 
				<< (obj.gap_indexs)[i].first
				<< " " 
				<< (obj.gap_indexs)[i].second 
				<< endl;
		}
	}

	LOG4CXX_DEBUG(logger, "Output each scaffold segments end...");
	return os;
}
