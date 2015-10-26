#include "gap_filler.h"
#include "component.h"
#include "condensed_debruijn_graph.h"
#include "condensed_debruijn_graph_reader.h"
#include "constant.h"
#include "contigs.h"
#include "kmer.h"

#include <deque>
#include <iterator>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
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
    LOG4CXX_DEBUG(logger, "Begin Gap Filling ... ");

	int cur_gap_index(0);

	for (size_t i = 0; i < _scaffolds.size(); ++i) {
		if (_scaffolds[i].contigs.size() <= 1) { // no gaps
			continue;
		}
        Component::ContigIdList contigs = _scaffolds[i].contigs;

		for (size_t j = 1; j < contigs.size(); ++j) {
			int left_index = contigs[j - 1];
			int right_index = contigs[j];

			gap_indexs.push_back(make_pair(left_index, right_index));

            const std::string& lsequence = _uniq_graph._indexer[left_index].seq;
            const std::string& rsequence = _uniq_graph._indexer[right_index].seq;
            std::string suffix = lsequence.substr(lsequence.length() - (K - 1), K - 1);
            std::string prefix = rsequence.substr(0, K - 1);

			int overlap = alignment(suffix, prefix);
			if (overlap >= OVERLAP) {
                _gapinfo_tbl[std::make_pair(i, j)] = GapInfo(-1, overlap);
			} else if (_scaffolds[i].gaps[j - 1] > MU + 3*var) {
                _gapinfo_tbl[std::make_pair(i, j)] = GapInfo();
			    //fail_gap_info.push_back(cur_gap_index);
            } else {
                int distance = _scaffolds[i].gaps[j - 1];
                try {
                    GapInfo gapinfo;
                    BFS(cur_gap_index, left_index, right_index, distance, &gapinfo);
                    _gapinfo_tbl[std::make_pair(i, j)] = gapinfo;
#ifdef RESULT 
                    //count_base_in_gap += distance;
#endif
                } catch(bad_alloc) {
                    LOG4CXX_WARN(logger, "BFS is too memory-intensive, ignoring...");
                    _gapinfo_tbl[std::make_pair(i, j)] = GapInfo();
                    //fail_gap_info.push_back(cur_gap_index);
                }
            }
			cur_gap_index++;
		}
	}
	LOG4CXX_DEBUG(logger, boost::format("The number of gaps = %d") % gap_indexs.size());
	//LOG4CXX_DEBUG(logger, boost::format("Filling the gaps which shares overlap ( < %d) = %d") % OVERLAP % pre_gap_info.size());

	assert(gap_indexs.size() == cur_gap_index);
	assert(_gapinfo_tbl.size() == cur_gap_index);
#ifdef  RESULT
    //LOG4CXX_DEBUG(logger, boost::format("Base counts in gaps : %d") % count_base_in_gap);
#endif
	//LOG4CXX_DEBUG(logger, boost::format("Unique candidate gaps number = %d") % uniq_gap_info.size());
	LOG4CXX_DEBUG(logger, boost::format("Multi candidates gaps number = %d") % multi_gap_info.size());
	//LOG4CXX_DEBUG(logger, boost::format("Failed gaps number = %d") % fail_gap_info.size());
	//LOG4CXX_DEBUG(logger, boost::format("Pre  gaps number = %d") % pre_gap_info.size());
	multi_candidates_gap_filling();
	output_initial_scaffolds_seq();
	get_multi_candidates_gap_local_seq();
}

//Align two nerghboring contigs
size_t GapFiller::alignment(const string& suffix, const string& prefix) {
	BOOST_ASSERT(prefix.length() == suffix.length());

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

struct _ContigLengthCmp_ {
    bool operator()(const std::string& c1, const std::string& c2) const {
        return c1.length() < c2.length();
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
        _scaffolds.push_back(component);
    }
    return true;
}

bool GapFiller::input_contigs(const std::string& file) {
	//get unique edges for training parameter
#ifdef TRAINING
    std::priority_queue< std::string, std::vector< std::string >, _ContigLengthCmp_ > unique_contigs_for_training;
#endif

    // cdbg_copy_number.fa
    {
        std::ifstream stream(file.c_str());
        if (!stream) {
            LOG4CXX_WARN(logger, boost::format("%s: No such file or directory!!!") % file);
            return false;
        }
        ContigReader reader(stream);
        Contig contig;
        while (reader.read(contig)) {
#ifdef TRAINING
            if (contig.copy_num == 1) {
                unique_contigs_for_training.push(contig.seq);
            }
#endif
            _uniq_graph.addEdge(contig.seq, contig.copy_num);
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

    return true;
}

std::string GapFiller::path2seq(const CondensedDeBruijnGraph& graph, const Path& path) const {
    std::string seq;
    if (!path.empty()) {
        seq += graph._indexer[path[0]].seq;
        for (size_t i = 1; i < path.size(); ++i) {
            seq += graph._indexer[path[i]].seq.substr(K - 1);
        }
    }
    return seq;
}

void GapFiller::BFS(const CondensedDeBruijnGraph& graph, const std::string& lseq, const std::string& rseq, int gap, size_t max_depth, size_t max_queue, PathList& pathlist) {

    std::string suffix = lseq.substr(lseq.length() - (K - 1), K - 1);
    std::string prefix = rseq.substr(0, K - 1);

    CondensedDeBruijnGraph::NodeList::const_iterator i = graph._children.find(suffix), j = graph._parents.find(prefix);
    if (i != graph._children.end() && j != graph._parents.end()) {
        for (CondensedDeBruijnGraph::EdgeList::const_iterator m = i->second.begin(); m != i->second.end(); ++m) {
            if (m->first == lseq[lseq.length() - K]) {
                for (CondensedDeBruijnGraph::EdgeList::const_iterator n = j->second.begin(); n != j->second.end(); ++n) {
                    if (n->first == rseq[K - 1]) {
                        BFS(graph, m->second, n->second, gap, STEP, 2000, pathlist);
                    }
                }
            }
        }
    }
}

void GapFiller::BFS(const CondensedDeBruijnGraph& graph, const size_t from, const size_t to, int gap, size_t max_depth, size_t max_queue, PathList& pathlist) {
	size_t step = 0;

    Path path = boost::assign::list_of(from);
    typedef std::pair< Path, int > Choice;
    std::deque< Choice > Q = boost::assign::list_of(std::make_pair(
                path, 0 - (K - 1)
                ));
    while (!Q.empty()) {
        ++step;
		if (step > max_depth) {
			LOG4CXX_DEBUG(logger, boost::format("BFS step is bigger than %d...") % max_depth);
			break;
		}

        if (Q.size() > max_queue) {
			LOG4CXX_DEBUG(logger, boost::format("BFS Q size=%d is greater than %d...") % Q.size() % max_queue);
            break;
        }

        std::deque< Choice > nextQ;
        while (!Q.empty()) {
            Choice choice = Q.front();
            Q.pop_front();

            size_t node = choice.first.back();
            if (node == to) {
                pathlist.push_back(choice.first);
            }
            if (choice.second > gap - (K - 1)) {
                continue;
            }
            const CondensedDeBruijnGraph::CondensedEdge& edge = graph._indexer[node];
            std::string suffix = edge.seq.substr(edge.seq.length() - (K - 1), K - 1);
            CondensedDeBruijnGraph::NodeList::const_iterator i = graph._parents.find(suffix);
            if (i != graph._parents.end()) {
                for (CondensedDeBruijnGraph::EdgeList::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
                    //size_t length = graph._indexer[j->second].seq.length();
                    size_t length = graph._indexer[j->second].seq.length() - (K - 1);

                    Choice new_choice(choice);
                    new_choice.first.push_back(j->second);
                    new_choice.second += length;

                    nextQ.push_back(new_choice);
                }
            }
        }

        Q = nextQ;
    }
}

// Run BFS to fill current gap 
void GapFiller::BFS(const int gap_index, const  int left_index, const  int right_index, int dis, GapInfo* gapinfo) {
	vector<vector<int> >  final_sets; //only contain the middle sequence except two ends
	dis += (K-1+3*var); //gap constraints
	dis += 30;

    {
        PathList pathlist;
        BFS(_uniq_graph, left_index, right_index, dis, STEP, 1000, pathlist);

        BOOST_FOREACH(const Path& path, pathlist) {
            std::vector< int > new_path;
            std::copy(path.begin(), path.end(), std::back_inserter(new_path));
            if (!new_path.empty()) new_path.erase(new_path.begin());
            if (!new_path.empty()) new_path.pop_back();
            final_sets.push_back(new_path);
        }
        if (!pathlist.empty()) {
            if (pathlist.size() <= MAX_CHOICE) {
                if (gapinfo != NULL) {
                    gapinfo->graph = 0;
                    gapinfo->pathlist = pathlist;
                }
            } else {
                if (gapinfo != NULL) {
                    gapinfo->graph = -1;
                }
            }
        }
    }

	if (!final_sets.empty()) {
		if (final_sets.size() == 1) {
			//uniq_gap_info.push_back(make_pair(gap_index,final_sets[0]));
		} else {
			assert(final_sets.size() > 1);
			if (final_sets.size() > MAX_CHOICE) {
				//fail_gap_info.push_back(gap_index);
			} else {
				assert(final_sets.size() <= MAX_CHOICE);
				multi_gap_info.push_back(make_pair(gap_index, final_sets));
			}
		}

		assert(gap_state.find(gap_index) == gap_state.end());
		gap_state.insert(gap_index);
		return;
	}

    {
        PathList pathlist;

        const std::string& lsequence = _uniq_graph._indexer[left_index].seq;
        const std::string& rsequence = _uniq_graph._indexer[right_index].seq;
        BFS(_all_graph, lsequence, rsequence, dis, STEP, 2000, pathlist);

        BOOST_FOREACH(const Path& path, pathlist) {
            std::vector< int > new_path;
            std::copy(path.begin(), path.end(), std::back_inserter(new_path));
            if (!new_path.empty()) new_path.erase(new_path.begin());
            if (!new_path.empty()) new_path.pop_back();
            final_sets.push_back(new_path);

        }
        if (!pathlist.empty()) {
            if (pathlist.size() <= MAX_CHOICE) {
                if (gapinfo != NULL) {
                    gapinfo->graph = 0;
                    gapinfo->pathlist = pathlist;
                }
            } else {
                if (gapinfo != NULL) {
                    gapinfo->graph = -1;
                }
            }
        }
    }

	if (!final_sets.empty())
	{
		if (final_sets.size() == 1)
		{
			//uniq_gap_info.push_back(make_pair(gap_index,final_sets[0]));
		}	
		else
		{
			if (final_sets.size() > MAX_CHOICE)
			{
				//fail_gap_info.push_back(gap_index);
			}
			else
			{
				assert(final_sets.size() <= MAX_CHOICE);
				multi_gap_info.push_back(make_pair(gap_index, final_sets));
			}
		}

		return;
	}

	//fail_gap_info.push_back(gap_index);
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

//get initial scaffolds with no filling gap filling and using 'N' to fill failed gaps
void GapFiller::output_initial_scaffolds_seq() {
	int cur_gap_index = 0 ;

	long GENOME = 0;		//estimate the genome size using total length of scaffold
    std::string file = boost::str(boost::format("%dmer.scaf_seq_with_gaps") % K);

    //(index , length) pair of scaffold
    std::vector< size_t > scaffold_length_list;
    ofstream out_scaf_seq(file.c_str());
    if (!out_scaf_seq) {
        cerr << "[Info] Create " << file << " failed..." << endl;
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < _scaffolds.size(); ++i) {
        std::string seq;         //scaffold sequence
        std::string one_segment; //one segment of scaffold since the multi gaps

        //for scaffolds segments
        std::vector< std::string > cur_scaff_segments;

        const Component::ContigIdList& contigs = _scaffolds[i].contigs;
        BOOST_ASSERT(!contigs.empty());
        seq += _uniq_graph._indexer[contigs[0]].seq;
        one_segment += _uniq_graph._indexer[contigs[0]].seq;

        if (contigs.size() == 1) {
            scaffold_length_list.push_back(seq.length());
            out_scaf_seq << ">scaf_" << i << "_" << seq.size()  << endl;
            out_scaf_seq << seq << endl;
            GENOME += seq.size();

            //for scaf segment
            cur_scaff_segments.push_back(one_segment);
            _scaffolds_segments.push_back(cur_scaff_segments);

            continue;
        }

        for (size_t j = 1; j < contigs.size(); ++j) {
            int left_index = contigs[j - 1];
            int right_index = contigs[j];

            int left_len = _uniq_graph._indexer[left_index].seq.length();

            GapInfoTable::iterator it = _gapinfo_tbl.find(std::make_pair(i, j));
            BOOST_ASSERT(it != _gapinfo_tbl.end());

            assert(gap_indexs[cur_gap_index].first == left_index && gap_indexs[cur_gap_index].second == right_index);

            if (it->second.graph == -1 && it->second.overlap > 0) {    // overlap
                BOOST_ASSERT(one_segment.size() > it->second.overlap);
                one_segment += _uniq_graph._indexer[right_index].seq.substr(it->second.overlap);
                BOOST_ASSERT(it->second.graph == -1);
            } else if (it->second.graph == -1) {                     // failed gap
                BOOST_ASSERT(it->second.graph == -1);
                BOOST_ASSERT(it->second.overlap == 0);
                //for failed gap , the count of 'N' is determined by the distance estimation
                int distance = _scaffolds[i].gaps[j - 1];
                if (distance > (MU + 3*var)) {
                    one_segment += std::string(MU+3*var,'N');
                } else if (distance < 0) {
                } else {
                    one_segment += std::string(distance, 'N');
                }
                one_segment += _uniq_graph._indexer[right_index].seq;
            } else if (it->second.graph != -1 && it->second.pathlist.size() == 1) { // unique
                assert(one_segment.size() >= left_len);
                one_segment.resize(one_segment.size() - left_len);
                one_segment += path2seq(it->second.graph == 0 ? _uniq_graph : _all_graph, it->second.pathlist[0]);
                BOOST_ASSERT(it->second.graph != -1);
            } else {                                                // multi_gap
                BOOST_ASSERT(it->second.pathlist.size() <= MAX_CHOICE);
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

        _scaffolds_segments.push_back(cur_scaff_segments);

        //initial scaff
        seq += one_segment;
        out_scaf_seq << ">scaf_" << i << "_" <<  seq.size()  << endl;
        out_scaf_seq << seq << endl;
        GENOME += seq.size();

        scaffold_length_list.push_back(seq.length());
    }
    assert(_scaffolds_segments.size() ==  _scaffolds.size());
    assert(_scaffolds.size() ==  scaffold_length_list.size());
    assert(initial_multi_gaps.size() == multi_gap_info.size());

    out_scaf_seq.close();
    out_scaf_seq.clear();

	//output result
    {
        std::string file = boost::str(boost::format("%dmer.scaf_len") % K);
        std::ofstream stream(file.c_str());
        if (!stream) {
            LOG4CXX_ERROR(logger, boost::format("%s No such file or directory!!!") % file);
            exit(EXIT_FAILURE);
        }

        //find top 20
        std::sort(scaffold_length_list.begin(), scaffold_length_list.end(), std::greater< size_t >());

        for (size_t i = 0; i < scaffold_length_list.size(); ++i) {
            stream << scaffold_length_list[i] << std::endl;
        }
    }

    LOG4CXX_INFO(logger, boost::format("Genome estimated size : %d") % GENOME);

	long N50 = GENOME * 0.5 ;
	long N90 = GENOME * 0.9 ;

	LOG4CXX_INFO(logger, "---Scaffold---");
	LOG4CXX_INFO(logger, "Scaff length : Top 20...");

	for (size_t i = 0; i < min((size_t)20, scaffold_length_list.size()); ++i) {
        LOG4CXX_INFO(logger, boost::format("%d") % scaffold_length_list[i]);
	}

    std::vector< size_t > sum (scaffold_length_list.size(), 0);
	
    if (!scaffold_length_list.empty()) {
        sum[0] = scaffold_length_list[0];
        for (size_t i = 1; i < scaffold_length_list.size(); ++i) {
            sum[i] +=  sum[i - 1] + scaffold_length_list[i];
        }
    }
	long N50_val = 0;
	long N90_val = 0;

	for (size_t i = 0 ; i < sum.size(); ++i) {
		if (sum[i] >= N50 && N50_val == 0) {
			N50_val = scaffold_length_list[i];
		}
		if (sum[i] >= N90 && N90_val == 0) {
			N90_val = scaffold_length_list[i];
			break;
		}
	}

	LOG4CXX_INFO(logger, boost::format("N50 : %d") % N50_val);
	LOG4CXX_INFO(logger, boost::format("N90 : %d") % N90_val);

	LOG4CXX_INFO(logger, boost::format("All gaps number in scaffolds : %d") % gap_indexs.size());
	//LOG4CXX_INFO(logger, boost::format("filled gap numbers : %d") % (gap_indexs.size() - fail_gap_info.size()));
	//LOG4CXX_INFO(logger, boost::format("failed gap numbers : %d") % fail_gap_info.size());
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
        for(int i=0;i<_scaffolds_segments.size() ;++i)
        {
            debug_count += (_scaffolds_segments[i].size()-1);
        }
        assert(debug_count == multi_gap_info.size());
    }

    //cout << _scaffolds_segments.size() << endl;

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

	for(int i = 0 ; i < _scaffolds_segments.size() ; ++i)
	{
		cur_scaffold_segment.clear();

		vector<string> seqs = _scaffolds_segments[i];

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
    /*
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
    */

	os << "each scaffold segment..." << std::endl;
	assert(obj._scaffolds_segments.size() == obj._scaffolds.size());
	for (size_t i = 0; i < obj._scaffolds_segments.size() ; ++i) {
		assert(obj._scaffolds_segments[i].size() > 0);

        os << boost::format("i = %d") % i << std::endl;

        std::copy(
                obj._scaffolds_segments[i].begin(), 
                obj._scaffolds_segments[i].end(),
                std::ostream_iterator< std::string >(os, "\n")
                );
		os << endl;
	}

	for (std::vector< std::vector< std::string > >::const_iterator 
			it = obj._scaffolds_segments.begin();
			it != obj._scaffolds_segments.end();
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
