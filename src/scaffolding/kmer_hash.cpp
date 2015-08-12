#include "kmer_hash.h"

#include <vector>
#include <fstream>
#include <algorithm>
#include <assert.h>

bool KmerHash::read(std::istream& stream) {
    return true;
}

bool KmerHash::read(const std::string& file) {
    std::ifstream stream(file.c_str());
    return read(stream);
}

using namespace std;

vector<Kmer> tmp_kmer_array;

void Kmer_Hash::push_a_kmer(Kmer tmp, long i, long j)
{
	if (kmer_node.kmer_num < 100000)
	{
		kmer_node.kmer_buff[kmer_node.kmer_num] = tmp;
		kmer_node.kmer_num ++;
		pos_node.pos_buff[pos_node.pos_num] = pair<int, int>(i, j);
		pos_node.pos_num ++;
	}else
	{
		kmer_buff.content.push_front(kmer_node);
		kmer_node.kmer_num = 0;
		kmer_node.kmer_buff[kmer_node.kmer_num] = tmp;
		kmer_node.kmer_num += 1;
		kmer_buff.kmer_num += 100000;

		pos_buff.content.push_front(pos_node);
		pos_node.pos_num = 0;
		pos_node.pos_buff[pos_node.pos_num] = pair<long, long>(i, j);
		pos_node.pos_num += 1;
		pos_buff.pos_num += 100000;
	}
}

void Kmer_Hash::free_buff()
{
	cout << "begin free kmer and position buffer" << endl;
	kmer_buff.content.clear();
	pos_buff.content.clear();
	tmp_kmer_array.clear();
}

unsigned long Kmer_Hash::size()
{
	return kmer_array.size();
}

void Kmer_Hash::clear()
{
	kmer_array.clear();
	pos_array.clear();
	kmer_index.clear();
}

bool mycmpkmer(long i, long j)
{
	if ((tmp_kmer_array[i] < tmp_kmer_array[j])
			||
			(tmp_kmer_array[i] == tmp_kmer_array[j] & tmp_kmer_array[i] < tmp_kmer_array[j])
			||
			(tmp_kmer_array[i] == tmp_kmer_array[j] & tmp_kmer_array[i] == tmp_kmer_array[j] & tmp_kmer_array[i] < tmp_kmer_array[j])
	   )
		return true;
	else
		return false;

};

void Kmer_Hash::initialize_kmer_array()
{
	cout << "begin initialize kmer array" << endl;
	unsigned long count = 0;
	kmer_array.resize(kmer_buff.kmer_num + kmer_node.kmer_num);
	pos_array.resize(kmer_buff.kmer_num + kmer_node.kmer_num);
	kmer_index.resize(kmer_buff.kmer_num + kmer_node.kmer_num);

	for(int i = 0; i < kmer_array.size(); i++)
	{
		kmer_index[i] = i;
	}

	list<Kmer_Node>::iterator it = kmer_buff.content.begin();
	cout << "\tkmer number = " << kmer_buff.kmer_num + kmer_node.kmer_num << endl;
	for (; it != kmer_buff.content.end() ; it++)
	{
		for (unsigned long j = 0; j < 100000 && j < it->kmer_num; j++)
		{
			kmer_array[count] = it->kmer_buff[j];
			count ++;
		}
	}
	for (long i = 0; i < kmer_node.kmer_num; i++)
	{
		kmer_array[count] = kmer_node.kmer_buff[i];
		count ++;
	}
	//cout << "the number of kmer copyed is " << count << endl;

	count = 0;

	list<Posi_Node>::iterator it_pos = pos_buff.content.begin();
	//cout << "pos number is " << pos_buff.pos_num + pos_node.pos_num << endl;
	for (; it_pos != pos_buff.content.end() ; it_pos++)
	{
		for (unsigned long j = 0; j < 100000 && j < it_pos->pos_num; j++)
		{
			pos_array[count] = it_pos->pos_buff[j];
			count ++;
		}
	}
	for (long i = 0; i < kmer_node.kmer_num; i++)
	{
		pos_array[count] = pos_node.pos_buff[i];
		count ++;
	}

	//cout << "the number of kmer copyed is " << count << endl;

	count = 0;
	tmp_kmer_array.resize(kmer_array.size());
	for(int i = 0; i < tmp_kmer_array.size(); i++)
	{
		tmp_kmer_array[i] = kmer_array[i];
		count ++;
	}

	//cout << "the number of kmer copyed is " << count << endl;

	//cout << "tmp kmer array size " << tmp_kmer_array.size() << endl;
	sort(kmer_index.begin(), kmer_index.end(), mycmpkmer);
	//cout << "sort kmer index end" << endl;



}

pair<long, long> Kmer_Hash::get_kmer_pos(Kmer tmp)
{
	long low, media, up;
	low = 0;
	up = kmer_array.size() - 1;
	media = (low + up) / 2;
	while (low <= up)
	{
		if (kmer_array[kmer_index[media]] == tmp)
		{
			return pos_array[kmer_index[media]];
		}else if (kmer_array[kmer_index[media]] < tmp)
		{
			low = media + 1;
		}else
		{
			up = media - 1;
		}
		media = (up + low) / 2;
	}
	return pair<long, long>(-1, -1);
}


