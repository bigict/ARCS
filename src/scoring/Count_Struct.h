#ifndef _COUNT_STRUCT_H
#define _COUNT_STRUCT_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "Constant.h"

using namespace std;

struct Count_Struct
{
	friend ostream &operator<<(ostream &, const Count_Struct &);
	
	vector< vector< long > > error_types;       //base mismatch
	vector< long > base_counts;                 //base number of A C T G N
	vector< long > read_lens;                   //count of read length
	vector< long > insert_lens;                 //count of insert length

	vector< int > error_pos;                    //error(indel, mismatch) pos in read
	vector< int > ins_pos;                      //count of insertion at each position of read
	vector< int > ins_len;                      //count of insertion length at each position of read
	vector< int > del_pos;                      //count of deletion at each position of read
	vector< int > del_len;                      //count of deletion length at each position of read

	vector< long > effective_lens;              //effective sequence lengths

	long unique_mapped_count;                   //unique mapped read count
	long unused_read_count;                     //mapped which is invalid
	long error_read_count;                      //error mapped read count
	long max_insert_len;                        //max insert size length
	int max_read_len;                           //max read length
	long unmapped_read_count;                   //not used here
	long total_read_count;                      //total read count
	
    //training contigs
	vector< string > contigs;                  
	vector< long > contig_lens;                 
	long total_contig_len ;

    //ctor
	Count_Struct():
		unique_mapped_count(0),
		unused_read_count(0),
		error_read_count(0),
		max_insert_len(0),
		max_read_len(0),
		unmapped_read_count(0),
		total_read_count(0),
		total_contig_len(0),
		error_types(ALPHABET, vector< long >(ALPHABET, 0)),
		base_counts(ALPHABET, 0)
	{
      //
	}
};

#endif  /*_COUNT_STRUCT_H*/
