#include "Kmer.h"

extern int K;

Kmer::Kmer()
{
}

Kmer::Kmer(string const &str)
{
	second_code = 0;
	for(int i = 0; i < 32; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 32; i < K; i++)
	{
		second_code = (second_code << 2) + ((str[i] & 0x06) >> 1);
	}

}

Kmer::Kmer(unsigned long first, unsigned long second)
{
	first_code = first;
	second_code = second;
}


Kmer Kmer::get_binary_code() const
{
	return *this;
}
unsigned long Kmer::get_code() const{
	return first_code | second_code;
}

void Kmer::set_seq(string const &str)
{
	second_code = 0;
	for(int i = 0; i < 32; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 32; i < K; i++)
	{
		second_code = (second_code << 2) + ((str[i] & 0x06) >> 1);
	}
}

extern string tem;

string Kmer::get_seq() const
{
	char ch[64];
	unsigned long tmp_first_code = first_code;
	unsigned long tmp_second_code = second_code;
	int tmp_i;
	for(int i = 0; i < 32; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[31-i] = tem[tmp_i];
	}
	for(int i = 32; i < K; i++)
	{
		tmp_i = tmp_second_code & 0x0000000000000003;
		tmp_second_code = tmp_second_code >> 2;
		ch[K - i + 31] = tem[tmp_i];
	}
	return string(ch, K);

}

Kmer_Short Kmer::get_next_node()
{
	unsigned long tmp_first_code = first_code;
	unsigned long tmp_second_code = second_code;
	tmp_first_code = tmp_first_code << 2;
	unsigned long ff = 0x0000000000000003;
	for (int i = 0; i < (K - 33); i++)
	{
		ff = ff << 2;
	}
	unsigned long last = ff & tmp_second_code;
	ff = ~ff;
	tmp_second_code &= ff; 
	for (int i = 0; i < (K - 33); i++)
	{
		last = last >> 2;
	}
	tmp_first_code += last;

	return Kmer_Short(tmp_first_code, tmp_second_code);
}

Kmer_Short Kmer::get_pre_node()
{
	unsigned long tmp_second_code = second_code;
	tmp_second_code = tmp_second_code >> 2;
	return Kmer_Short(first_code, tmp_second_code); 
}

Kmer &Kmer::operator = (const string & str)
{
	second_code = 0;
	for(int i = 0; i < 32; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 32; i < K; i++)
	{
		second_code = (second_code << 2) + ((str[i] & 0x06) >> 1);
	}
	return *this;
}

bool Kmer::operator == (const Kmer &tmp)
{
	if(first_code == tmp.first_code & second_code == tmp.second_code)
		return true;
	else
		return false;
}

bool Kmer::operator != (const Kmer &tmp)
{
	if(first_code != tmp.first_code || second_code != tmp.second_code)
		return true;
	else
		return false;
}

bool Kmer::operator < (const Kmer &tmp)
{
	if(first_code < tmp.first_code || (first_code == tmp.first_code & second_code < tmp.second_code))
		return true;
	else
		return false;
}

bool Kmer::operator > (const Kmer &tmp)
{
	if(first_code > tmp.first_code || (first_code == tmp.first_code & second_code > tmp.second_code))
		return true;
	else
		return false;
}

ostream &operator << (ostream &out, const Kmer &tmp)
{
	char ch[64];
	unsigned long tmp_first_code = tmp.first_code;
	unsigned long tmp_second_code = tmp.second_code;
	int tmp_i;
	for(int i = 0; i < 32; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[31 - i] = tem[tmp_i];
	}
	for(int i = 32; i < K; i++)
	{
		tmp_i = tmp_second_code & 0x0000000000000003;
		tmp_second_code = tmp_second_code >> 2;
		ch[K - i + 31] = tem[tmp_i];
	}
	out << string(ch, K);
	return out;
}
