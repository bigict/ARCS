#include "Para.h"

ostream &operator<<(ostream &os, const Para &para)
{
	os << "Unique mapped read count\t=\t" << para.unique_mapped_count << endl;
	os << "unused read count\t=\t: " << para.unused_read_count << endl;
	os << "error read count\t=\t" << para.error_read_count << endl;
	os << "max insert len\t=\t" << para.max_insert_len << endl;
	os << "max read len\t=\t: " << para.max_read_len << endl;

	os << "A\tT\tC\tG\tN" << endl;

	for(int i = 0 ; i < ALPHABET ;++i)
	{
		for(int j = 0 ; j <  ALPHABET; ++j)
		{
			os << para.base_error_type[i][j] << "\t";
		}
		os << endl;
	}
	
	os << endl;

	os << "base_error_rate : "<< endl;
	for(int i = 0 ; i < ALPHABET ;++i)
		os << para.base_error_rate[i] << "\t";
	os << endl;

	os << endl;

	int index = 0;

	while(index < para.max_read_len)
	{
		os << "err_pos[" << index << "] = " << para.error_pos_dist[index] << "\t";
		os << "ins_pos[" << index << "] = " << para.ins_pos_dist[index] << "\t";
		os << "ins_len[" << index << "] = " << para.ins_len_dist[index] << "\t";
		os << "del_pos[" << index << "] = " << para.del_pos_dist[index] << "\t";
		os << "del_len[" << index << "] = " << para.del_len_dist[index] << "\t";
		os << "no_err[" << index << "] = " << para.no_error_prob[index] << "\t";
		os << endl;
		index++;
	}
	os << endl << endl;



	index = 0;
	while(index < para.max_insert_len)
	{
		os << index << "\t" << para.insert_len_dist[index] << "\t";
		os << "effective_lens[" << index << "]\t=\t" << para.effective_lens[index] << endl;
		index++;
	}

	os << endl << endl;

	os << "total_contig_len\t=\t" << para.total_contig_len << endl;
	for(int i = 0 ; i < para.contig_lens.size() ; ++i)
	{
		os << para.contig_lens[i] << endl;
	}
	os << endl;

	//copy(para.contig_lens.begin(), para.contig_lens.end(), ostream_iterator<long>(os , "\n"));
    //os << endl;
	return os;
}
