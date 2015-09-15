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
#include "debruijn_edge.h"

using namespace std;

extern int K;
extern int overlap;
extern int MU;
extern int var;
extern int EXTEND;
extern int STEP;

class Gap_Filling {
private:
	friend std::ostream& operator<<(std::ostream& os, const Gap_Filling& obj);

	typedef std::map< std::string, std::vector< int> > candidates_type;
	typedef std::vector< std::pair< std::vector< int >, int > > choice_type;
    typedef std::vector< Edge > EdgeList;


    EdgeList all_cdbg_edges;    //condensed contigs
    EdgeList all_initial_edges; //condensed contigs in initial CDBG
    ComponentList all_scaffolds; //scaffolds
    std::vector< int > gap_distances;      //gap distances in scaffolds

    //record the (gap, index) pair
    std::map< std::pair< int, int > , int > gap_indexs_hash;
    std::vector< std::pair< int, int > > gap_indexs;

	//record the gap is filled by condensed contigs or initial contigs 
    //true -- condensed contigs
	//false -- initial contigs
    std::set< int > gap_state;

    //temp data to reconstruct the CDBG
	candidates_type heads_initial_hash;
	candidates_type tails_initial_hash;
	
	candidates_type heads_cdbg_hash;
	candidates_type tails_cdbg_hash;

    //construct the CDBG
    std::vector< std::vector< int > > next_candidates_in_initial;
    std::vector< std::vector< int > > prev_candidates_in_initial;

	vector<vector<int> > next_candidates_in_cdbg;
	vector<vector<int> > prev_candidates_in_cdbg;

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

public:
	Gap_Filling() {
    }
	virtual ~Gap_Filling() {
    }

	void gap_filling();

    bool input_scaffold(const std::string& file);
    bool input_contigs(const std::string& file);
    bool input_debruijn(const std::string& file);

private:
	void get_gap_info();                                //get gap info

	void uniq_candidate_gap_filling();                  //fill unique gaps
	void multi_candidates_gap_filling();                //fill multi gaps

	void get_multi_candidates_gap_local_seq();          //get multi gaps local sequences

    //set function
    void build_condensed_de_bruijn_graph();

	int alignment(const string &, const string &);
    void merge_segments(const set<int > &);

	void get_all_scaffolds_segments();
    void get_gap_candidate_number_hash();

	void output_initial_scaffolds_seq();   // no multi gaps and filled 'N' for fail gap
	void BFS(const  int, const  int , const  int , int);
	void DFS(vector<int> &, const vector<int> &, vector<string> &, int , vector<string> &, string , string& , vector<string> &);
};

#endif /*_GAP_FILLING_H*/
