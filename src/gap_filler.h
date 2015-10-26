#ifndef gap_filler_h_
#define gap_filler_h_

#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <string>
#include <map>
#include <set>
#include <iterator>
#include <cstdlib>
#include <algorithm>
#include <queue>


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "component.h"
#include "condensed_debruijn_graph.h"

using namespace std;

extern int K;
extern int OVERLAP;
extern int MU;
extern int var;
extern int EXTEND;
extern int STEP;

class GapFiller {
public:
	GapFiller() : _uniq_graph(K), _all_graph(K) {
    }
	virtual ~GapFiller() {
    }

	void fill();

    bool input_scaffold(const std::string& file);
    bool input_contigs(const std::string& file);
    bool input_debruijn(const std::string& file);

private:
	friend std::ostream& operator<<(std::ostream& os, const GapFiller& obj);

	void multi_candidates_gap_filling();                //fill multi gaps

	void get_multi_candidates_gap_local_seq();          //get multi gaps local sequences
	void output_initial_scaffolds_seq();   // no multi gaps and filled 'N' for fail gap

	size_t alignment(const string& suffix, const string& prefix);

    typedef std::vector< size_t > Path;
    typedef std::vector< std::vector< size_t > > PathList;
    struct GapInfo {
        GapInfo(const PathList& pathlist, size_t graph, int overlap) : pathlist(pathlist), graph(graph), overlap(overlap) {
        }
        GapInfo(size_t graph, int overlap) : graph(graph), overlap(overlap) {
        }
        GapInfo() : graph(-1), overlap(0) {
        }
        PathList pathlist;
        size_t graph;     // graph index
        int overlap;      // overlap
    };
    typedef std::pair< size_t, size_t > GapIndex;
    typedef std::map< GapIndex, GapInfo > GapInfoTable;

    std::string path2seq(const CondensedDeBruijnGraph& graph, const Path& path) const;
    void BFS(const CondensedDeBruijnGraph& graph, const size_t i, const size_t j, int distance, size_t max_depth, size_t max_queue, PathList& pathlist);
    void BFS(const CondensedDeBruijnGraph& graph, const std::string& lseq, const std::string& rseq, int distance, size_t max_depth, size_t max_queue, PathList& pathlist);
	void BFS(const  int, const  int , const  int , int, GapInfo* gapinfo);

	void DFS(vector<int> &, const vector<int> &, vector<string> &, int , vector<string> &, string , string& , vector<string> &);

    //record the (gap, index) pair
    std::vector< std::pair< int, int > > gap_indexs;

	//record the gap is filled by condensed contigs or initial contigs 
    //true -- condensed contigs
	//false -- initial contigs
    std::set< int > gap_state;

	map< int, vector< string > > multi_gap_seq;                 //candidate seq (do not include left and right contigd)

    //gap classification
	vector< pair<int,vector<vector<int> > > > multi_gap_info;   // multi candidates
	

    //the remaind multi gaps
	vector<int> initial_multi_gaps;
    
    //the remaind multi gaps and index them again, 
    //the index is the initial index!!!
	map<int, vector< string > > multi_gap_candidate_seq;

	//initial scaffolds seq and several segs since the multi gap
	//include the left and right (K-1)-mer
    //Do not merge the gaps whose choices is grt MAX_CHOICE
	vector<vector<string> > _scaffolds_segments;

    ComponentList _scaffolds; //scaffolds
    CondensedDeBruijnGraph _uniq_graph;
    CondensedDeBruijnGraph _all_graph;

    GapInfoTable _gapinfo_tbl;
};

#endif //gap_filler_h_
