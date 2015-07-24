
#include "Count_Struct.h"

ostream &operator<<(ostream &os, const Count_Struct &cs)
{
    
	os << "Max read len\t=\t" << cs.max_read_len << endl;
	os << "Read lenths\t:\t" << endl;
	for(int i = 0 ; i < cs.read_lens.size() ; ++i)
	{
		os << i << "\t" << cs.read_lens[i] << endl;
	}

	os << "Effective segment lengths\t:\t" << endl;
	for(int i = 0 ; i < cs.effective_lens.size(); ++i)
	{
		os << i << "\t" << cs.effective_lens[i] << endl;
	}

	os << "Unmapped read count\t=\t" << cs.unmapped_read_count << endl;
	os << "Total read count\t=\t" << cs.total_read_count << endl;
	os << "Unique mapped read count\t=\t" << cs.unique_mapped_count << endl;
	os << "Unused read count\t=\t" << cs.unused_read_count << endl;
	os << "Error read count\t=\t" << cs.error_read_count << endl;
	os << "Max insert size\t=\t" << cs.max_insert_len << endl;

	os << "Insert size count" << endl;
	for(int i = 0 ; i < cs.insert_lens.size(); ++i)
	{
		os << i << "\t" << cs.insert_lens[i] << endl;
	}

	os << "Error_types" << endl;
	for(int i = 0 ;i < ALPHABET;++i)
	{
		for(int j = 0 ; j < ALPHABET;++j)
		{
			os << BASE[i] << "\t" << BASE[j] << "\t" << cs.error_types[i][j] << endl; 
		}
	}
	os << "Base_counts" << endl;
	for(int i = 0 ;i < ALPHABET;++i)
	{
		os << BASE[i] << "\t" << cs.base_counts[i] << endl;
	}
	int index = 0;

	while(index < cs.max_read_len)
	{
		os << "err_pos[" << index << "]\t=\t" << cs.error_pos[index] << "\t";
		os << "ins_pos[" << index << "]\t=\t" << cs.ins_pos[index] << "\t";
		os << "ins_len[" << index << "]\t=\t" << cs.ins_len[index] << "\t";
		os << "del_pos[" << index << "]\t=\t" << cs.del_pos[index] << "\t";
		os << "del_len[" << index << "]\t=\t" << cs.del_len[index] << "\t";
		os << endl;
		index++;
	}
	os << endl;
	os << "Total contig length\t=\t" << cs.total_contig_len << endl;
	return os;
}
