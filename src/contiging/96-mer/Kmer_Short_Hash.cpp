#include "Kmer_Short_Hash.h"

#include <vector>
#include <iostream>
//#include <fstream>
#include <algorithm>

void Kmer_Short_Hash::push_a_kmer_short(Kmer_Short tmp)
{
	if ( short_buff_node.kmer_short_num < 100000)
	{
		short_buff_node.kmer_short_buff[short_buff_node.kmer_short_num] = tmp;
		short_buff_node.kmer_short_num ++;

	}else
	{
//		cout << "build kmer buffer num: " << kmer_short_buff.content.size() + 1 << endl;
		kmer_short_buff.content.push_front(short_buff_node);
		short_buff_node.kmer_short_num = 0;
		short_buff_node.kmer_short_buff[short_buff_node.kmer_short_num] = tmp;
		short_buff_node.kmer_short_num += 1;
		kmer_short_buff.kmer_short_num += 100000;

	}
}

bool mycmpkmershort(Kmer_Short tmp1, Kmer_Short tmp2)
{
	if ((tmp1.first_code < tmp2.first_code)
			||
			(tmp1.first_code == tmp2.first_code & tmp1.second_code < tmp2.second_code)
			||
			(tmp1.first_code == tmp2.first_code & tmp1.second_code == tmp2.second_code & tmp1.third_code < tmp2.third_code))
		return true;
	else
		return false;

}

void Kmer_Short_Hash::free_buff()
{
	cout << "Release kmer buffer" << endl;
	kmer_short_buff.content.clear();
}

unsigned long Kmer_Short_Hash::size()
{
	return kmer_short_array.size();
}

void Kmer_Short_Hash::clear()
{
	kmer_short_array.clear();
}

void Kmer_Short_Hash::initialize_kmer_short_array()
{
	cout << "Building (k-1)mer array" << endl;
	vector<Kmer_Short> tmp_kmer_short_array;
	tmp_kmer_short_array.resize(kmer_short_buff.kmer_short_num + short_buff_node.kmer_short_num);
	unsigned int count = 0;
	list<Short_Buff_Node>::iterator it = kmer_short_buff.content.begin();
	//cout << "list node number is " << kmer_short_buff.content.size() << endl;
	//cout << "kmer short buff kmer short num is " << kmer_short_buff.kmer_short_num << endl;
	//cout << "short buff node kmer short num is " << short_buff_node.kmer_short_num << endl;
	//cout << "\t(k-1)mer number = " << kmer_short_buff.kmer_short_num + short_buff_node.kmer_short_num  << endl;
	for (; it != kmer_short_buff.content.end() ; it++)
	{
		for (unsigned int j = 0; j < 100000 && j < it->kmer_short_num; j++)
		{
			tmp_kmer_short_array[count] = it->kmer_short_buff[j];
			count ++;
		}
	}
	
	for (int i = 0; i < short_buff_node.kmer_short_num; i++)
	{
		tmp_kmer_short_array[count] = short_buff_node.kmer_short_buff[i];
		count ++;
	}
	
	free_buff();
	//cout << "the number of kmer short copyed is " << count << endl;
	//cout << "begin sort" << endl;
	sort(tmp_kmer_short_array.begin(), tmp_kmer_short_array.end(), mycmpkmershort);
	//cout << "sort end" << endl;
	unsigned int kmer_short_count = 1;
	Kmer_Short cur = tmp_kmer_short_array[0];
	for (int i = 1; i < tmp_kmer_short_array.size(); i ++)
	{
		if (tmp_kmer_short_array[i] != cur)
		{
			cur = tmp_kmer_short_array[i];
			kmer_short_count ++;
		}
	}
	
	cout << "\t(k-1)mer number = " << kmer_short_count << endl;

	kmer_short_array.resize(kmer_short_count);

	unsigned int kmer_short_index = 1;
	unsigned int pos = 0;
	Kmer_Short cur_kmer = tmp_kmer_short_array[pos];
	kmer_short_array[0] = tmp_kmer_short_array[0];
	while (pos < tmp_kmer_short_array.size())
	{
		while (cur_kmer == tmp_kmer_short_array[pos] &	pos < tmp_kmer_short_array.size())
		{
			pos ++;
		}
		
		if (pos >= tmp_kmer_short_array.size())
		{
			break;
		}

		kmer_short_array[kmer_short_index] = tmp_kmer_short_array[pos];
		cur_kmer = tmp_kmer_short_array[pos];
		kmer_short_index ++;
	}


	tmp_kmer_short_array.clear();
}

long Kmer_Short_Hash::get_kmer_short_pos(Kmer_Short tmp)
{
	long low, media, up;
	low = 0;
	up = kmer_short_array.size() - 1;
	media = (low + up) / 2;
	while (low <= up)
	{
		if (kmer_short_array[media] == tmp)
		{
			return media;
		}else if (kmer_short_array[media] < tmp)
		{
			low = media + 1;
		}else
		{
			up = media - 1;
		}
		media = (up + low) / 2;
	}
	return -1;
}

Kmer_Short& Kmer_Short_Hash::operator[] (unsigned long pos)
{
	return kmer_short_array[pos];
}

Kmer_Short& Kmer_Short_Hash::at(unsigned long pos)
{
	return kmer_short_array.at(pos);
}

