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
//#include <unordered_set>


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Edge.h"
#include "Scaffold.h"
#include "Constant.h"

using namespace std;

extern int K;
extern int overlap;
extern int MU;
extern int var;
extern int EXTEND;
extern int STEP;
extern string EDGE_FILE_NAME;
extern string SCAFFOLD_FILE_NAME;
extern string TEMP_FILE_DIR;
extern string INITIAL_EDGE_FILE_NAME;

class Gap_Filling
{

	friend ostream &operator<<(ostream &, const Gap_Filling & );

private:

	typedef map< string, vector< int> > candidates_type;
	typedef map< string, vector< int> >::iterator candidates_iterator;
	typedef map< string, vector< int> >::value_type candidates_value_type;
	typedef map< string, vector< int> >::size_type candidates_size_type;

	typedef vector< pair< vector< int>, int> > choice_type;
	typedef vector< pair< vector< int>, int> >::iterator choice_value_iterator;
	typedef vector< pair< vector< int>, int> >::value_type choice_value_type;
	typedef vector< pair< vector< int>, int> >::size_type choice_size_type;


    //Input file name 
	string scaffold_file_name;        //scaffold file name
	string edge_file_name;            //contig file name
	string initial_edge_file_name;    //condensed contig file name


	vector< Edge > all_cdbg_edges;    //condensed contigs
	vector< Edge > all_initial_edges; //condensed contigs in initial CDBG
	vector< Scaffold > all_scaffolds; //scaffolds
    vector< int > gap_distances;      //gap distances in scaffolds

    //record the (gap, index) pair
	map<pair<int, int> , int> gap_indexs_hash;
	vector<pair<int,int> > gap_indexs;

	//record the gap is filled by condensed contigs or initial contigs 
    //true -- condensed contigs
	//false -- initial contigs
	set<int> gap_state;

    //temp data to reconstruct the CDBG
	candidates_type heads_initial_hash;
	candidates_type tails_initial_hash;
	
	candidates_type heads_cdbg_hash;
	candidates_type tails_cdbg_hash;

    //construct the CDBG
	vector<vector<int> > next_candidates_in_initial;
	vector<vector<int> > prev_candidates_in_initial;

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
	vector<pair<int,int> > new_multi_gaps;
	map<pair<int,int>,int> new_multi_gaps_hash;
    
    //the remaind multi gaps and index them again, 
    //the index is the initial index!!!
	map<int, vector< string > > multi_gap_candidate_seq;
	//map<pair<int,int> , vector<string> > new_multi_gap_candidate_seq;

	//initial scaffolds seq and several segs since the multi gap
	//include the left and right (K-1)-mer
    //Do not merge the gaps whose choices is grt MAX_CHOICE
	vector<vector<string> > all_scaffolds_segments;

	//temp data to scoring
	vector< vector< string> > new_all_scaffolds_segments;

public:

	Gap_Filling(const string & ,const string &, const string &);
	~Gap_Filling();

	void gap_filling();

private:
	void get_gap_info();                                //get gap info

	void uniq_candidate_gap_filling();                  //fill unique gaps
	void multi_candidates_gap_filling();                //fill multi gaps

	void get_multi_candidates_gap_local_seq();          //get multi gaps local sequences

    //set function
	void set_edges();
	void set_scaffolds();
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
