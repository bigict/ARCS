#include "Kmer_Hash.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

void Kmer_Hash::push_a_kmer(Kmer tmp)
{
	if (buff_node.kmer_num < 100000)
	{
		buff_node.kmer_buff[buff_node.kmer_num] = tmp;
		buff_node.kmer_num ++;
	}else
	{
		
		kmer_buff.content.push_front(buff_node);
		buff_node.kmer_num = 0;
		buff_node.kmer_buff[buff_node.kmer_num] = tmp;
		buff_node.kmer_num += 1;
		kmer_buff.kmer_num += 100000;
	}
}

void Kmer_Hash::free_buff()
{
	//cout << "begin destroy kmer buff" << endl;
	kmer_buff.content.clear();
}

unsigned long Kmer_Hash::size()
{
	return kmer_array.size();
}

void Kmer_Hash::clear()
{
	kmer_array.clear();
	kmer_num.clear();
}

bool mycmpkmer(Kmer tmp1, Kmer tmp2)
{
	if (tmp1.first_code < tmp2.first_code)
		return true;
	else if(tmp1.first_code == tmp2.first_code)
		return tmp1.second_code < tmp2.second_code;
	else
		return false;
}

void Kmer_Hash::initialize_kmer_array()
{
	cout << "Initializing kmer hash table" << endl;
	vector<Kmer> tmp_kmer_array;
	tmp_kmer_array.resize(kmer_buff.kmer_num + buff_node.kmer_num);
	unsigned int count = 0;
	list<Buff_Node>::iterator it = kmer_buff.content.begin();
	//cout << "kmer number is " << kmer_buff.kmer_num + buff_node.kmer_num << endl;
	for (; it != kmer_buff.content.end() ; it++)
	{
		for (unsigned int j = 0; j < 100000 && j < it->kmer_num; j++)
		{
			tmp_kmer_array[count] = it->kmer_buff[j];
			count ++;
		}
	}
	for (int i = 0; i < buff_node.kmer_num; i++)
	{
		tmp_kmer_array[count] = buff_node.kmer_buff[i];
		count ++;
	}
	
	free_buff();
//	cout << "the number of kmer copyed is " << count << endl;
//	cout << "begin sort kmer" << endl;
	
	sort(tmp_kmer_array.begin(), tmp_kmer_array.end(), mycmpkmer);
///	cout << "sort kmer end" << endl;

	unsigned long kmer_count = 0;	
	unsigned int cov = 0;
	unsigned long pos = 0;
	Kmer cur_kmer = tmp_kmer_array[pos];

	while (pos < tmp_kmer_array.size())
	{
		cov = 1;
		while (cur_kmer == tmp_kmer_array[pos + cov] && pos + cov < tmp_kmer_array.size())
		{
			cov ++;
		}
		if (cov > 1)
		{
			kmer_count ++;
		}
		cur_kmer = tmp_kmer_array[pos + cov];
		pos = pos + cov;
	}

	cout << "\tk-mer number = " << kmer_count << endl;

	kmer_array.resize(kmer_count);
	kmer_num.resize(kmer_count);
	
	//cout << "kmer array kmer num initialize end" << endl;

	unsigned long kmer_index = 0;
	pos = 0;
	cur_kmer = tmp_kmer_array[pos];

	while (pos < tmp_kmer_array.size())
	{
		cov = 1;
		while (cur_kmer == tmp_kmer_array[pos + cov] && pos + cov < tmp_kmer_array.size())
		{
			cov ++;
		}
		if (cov > 1)
		{
			kmer_array[kmer_index] = tmp_kmer_array[pos];
			kmer_num[kmer_index] = cov;	
			kmer_index ++;
		}
		cur_kmer = tmp_kmer_array[pos + cov];
		pos = pos + cov;
	}

	//cout << "kmer uniq count is " << kmer_count << endl;

	tmp_kmer_array.clear();
}

unsigned int Kmer_Hash::get_kmer_num(Kmer tmp)
{
	long low, media, up;
	low = 0;
	up = kmer_array.size() - 1;
	media = (low + up) / 2;
	while (low <= up)
	{
		if (kmer_array[media] == tmp)
		{
			return kmer_num[media];
		}else if (kmer_array[media] < tmp)
		{
			low = media + 1;
		}else
		{
			up = media - 1;
		}
		media = (up + low) / 2;
	}
	return 0;
}

unsigned int Kmer_Hash::get_kmer_num( unsigned long pos)
{
	return kmer_num[pos];
}

Kmer& Kmer_Hash::at(unsigned long pos)
{
	return kmer_array.at(pos);
}

Kmer& Kmer_Hash::operator [](unsigned long pos)
{
	return kmer_array[pos];
}

void Kmer_Hash::filter()
{
	cout << "Filtering out k-mers with coverage < 2" << endl;
	vector<Kmer> tmp_kmer_array = kmer_array;
	vector<unsigned int> tmp_kmer_num = kmer_num;
	
	unsigned long sum = 0;
	for(unsigned long i = 0; i < tmp_kmer_num.size(); i++)
	{
		if (tmp_kmer_num[i]  > 1)
		{
			sum ++;
		}
	} 
	
	kmer_array.resize(sum);
	kmer_num.resize(sum);
	
	sum = 0;
	for(unsigned long i = 0; i < tmp_kmer_num.size(); i ++)
	{
		if (tmp_kmer_num[i] > 1)
		{
			kmer_array[sum] = tmp_kmer_array[i];
			kmer_num[sum] = tmp_kmer_num[i];
			sum ++;
		}
	}
	
	cout << "\tk-mer number = " << sum << endl;
	tmp_kmer_array.clear();
	tmp_kmer_num.clear();
		
}
