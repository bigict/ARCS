#ifndef _GAP_FILLING_H
#define _GAP_FILLING_H

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
extern int overlap;
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

	void get_gap_info();                                //get gap info

	void uniq_candidate_gap_filling();                  //fill unique gaps
	void multi_candidates_gap_filling();                //fill multi gaps

	void get_multi_candidates_gap_local_seq();          //get multi gaps local sequences

	size_t alignment(const string& suffix, const string& prefix);
    void merge_segments(const set<int > &);

	void get_all_scaffolds_segments();
    void get_gap_candidate_number_hash();

    typedef std::vector< size_t > Path;
    typedef std::vector< std::vector< size_t > > PathList;

	void output_initial_scaffolds_seq();   // no multi gaps and filled 'N' for fail gap
    void BFS(const Component& component, std::vector< std::vector< size_t > >& pathlist);
    void BFS(const size_t i, const size_t j, std::vector< std::vector< size_t > >& pathlist);
	void BFS(const  int, const  int , const  int , int);
	void DFS(vector<int> &, const vector<int> &, vector<string> &, int , vector<string> &, string , string& , vector<string> &);

    ComponentList all_scaffolds; //scaffolds
    std::vector< int > gap_distances;      //gap distances in scaffolds

    //record the (gap, index) pair
    std::vector< std::pair< int, int > > gap_indexs;

	//record the gap is filled by condensed contigs or initial contigs 
    //true -- condensed contigs
	//false -- initial contigs
    std::set< int > gap_state;

    //gap (contig pair) -- gap candidate 
    map<pair<int,int> ,int > gap_candidate_number_hash;         //0-overlap, 1-uniq, >=2 :multi, 0-failed

	map< int, string > uniq_gap_seq;                            //candidate seq (including the left and right contigs)
	map< int, vector< string > > multi_gap_seq;                 //candidate seq (do not include left and right contigd)

    //gap classification
	vector< pair<int,vector<int> > > uniq_gap_info;             // unique candidate
	vector< pair<int,vector<vector<int> > > > multi_gap_info;   // multi candidates
	map< int , int> pre_gap_info;                               //overlap
	vector< int > fail_gap_info;                                //no candidate
	

    //the remaind multi gaps
	vector<int> initial_multi_gaps;
    
    //the remaind multi gaps and index them again, 
    //the index is the initial index!!!
	map<int, vector< string > > multi_gap_candidate_seq;
	//map<pair<int,int> , vector<string> > new_multi_gap_candidate_seq;

	//initial scaffolds seq and several segs since the multi gap
	//include the left and right (K-1)-mer
    //Do not merge the gaps whose choices is grt MAX_CHOICE
	vector<vector<string> > all_scaffolds_segments;

    CondensedDeBruijnGraph _uniq_graph;
    CondensedDeBruijnGraph  _all_graph;
};

#endif /*_GAP_FILLING_H*/
