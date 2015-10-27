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

    typedef std::vector< size_t > Path;
    typedef std::vector< std::vector< size_t > > PathList;
    struct GapInfo {
        GapInfo(const PathList& pathlist, size_t graph, int gap) : pathlist(pathlist), graph(graph), gap(gap) {
        }
        GapInfo(size_t graph, int gap) : graph(graph), gap(gap) {
        }
        GapInfo() : graph(-1), gap(0) {
        }
        PathList pathlist;
        size_t graph;       // graph index
        int gap;            // gap 
    };
    typedef std::pair< size_t, size_t > GapIndex;
    typedef std::map< GapIndex, GapInfo > GapInfoTable;

	size_t alignment(const string& suffix, const string& prefix);
    std::string path2seq(const CondensedDeBruijnGraph& graph, const Path& path) const;
    std::string path2seq(const CondensedDeBruijnGraph& graph, const Path& path, size_t i, size_t j) const;
    void BFS(const CondensedDeBruijnGraph& graph, const size_t i, const size_t j, int gap, size_t max_depth, size_t max_queue, PathList& pathlist);
    void BFS(const CondensedDeBruijnGraph& graph, const std::string& lseq, const std::string& rseq, int gap, size_t max_depth, size_t max_queue, PathList& pathlist);
	void BFS(size_t left_index, size_t right_index, int gap, GapInfo* gapinfo);

    ComponentList _scaffolds; //scaffolds
    CondensedDeBruijnGraph _uniq_graph;
    CondensedDeBruijnGraph _all_graph;

    GapInfoTable _gapinfo_tbl;
};

#endif //gap_filler_h_
