#ifndef _PARA_H
#define _PARA_H

#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <string>

#include <iterator>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>


#include "Constant.h"
#include "Count_Struct.h"

using namespace std;

struct Para
{
	
	friend ostream &operator<<(ostream &, const Para &);
	
	vector<vector<double> > base_error_type;      //error(mismatch) rate
	vector<double> base_error_rate;               //base error rate

	vector<double> error_pos_dist;                //error position distribution
	vector<double> ins_pos_dist;                  //insertion distribution 
	vector<double> ins_len_dist;                  //insertion length distribution
	vector<double> del_pos_dist;                  //deletion position distribution
	vector<double> del_len_dist;                  //deletion length distribution
	vector<double> insert_len_dist;               //insert size distribution
	vector<double> no_error_prob;                 //no error distribution
	vector< long > effective_lens;                //effecrive segments length distribution

	long unique_mapped_count;                     //unique mapped read count
	long unused_read_count;                       //unused read count
	long error_read_count;                        //error read count
	long max_insert_len;                          //max insert size
	int max_read_len;                             //max read length

	//scoring contigs
	vector<string> contigs;
	vector<long> contig_lens;
	long total_contig_len;

    //ctor
	Para():
		unique_mapped_count(0),
		unused_read_count(0),
		error_read_count(0),
		max_insert_len(0),
		max_read_len(0),
		base_error_type(ALPHABET, vector<double>(ALPHABET, 0)),
		base_error_rate(ALPHABET, 0),
		total_contig_len(0)
	{
      //
	}
};


#endif/*_PARA_H*/
