#include "Kmer_Short.h"

Kmer_Short::Kmer_Short()
{
}

Kmer_Short::Kmer_Short(string const &str)
{
	second_code = 0;
	for(int i = 0; i < 32; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 32; i < K - 1; i++)
	{
		second_code = (second_code << 2) + ((str[i] & 0x06) >> 1);
	}
}

Kmer_Short::Kmer_Short(unsigned long first, unsigned long second)
{
	first_code = first;
	second_code = second;
}
unsigned long Kmer_Short::get_code() const{
	return first_code | second_code;
}

Kmer Kmer_Short::get_next_kmer(unsigned int i) const
{
	unsigned long first = first_code;
	unsigned long second = second_code;
	second = second << 2;
	second += i;
	return Kmer(first, second);
}

Kmer_Short Kmer_Short::get_next_short(unsigned int i) const
{
	unsigned long second = second_code << 2;
	second |= i;
	return Kmer(first_code, second).get_next_node();
}

Kmer_Short Kmer_Short::get_next_short() const
{
	return Kmer(first_code, second_code).get_next_node();
}

Kmer_Short Kmer_Short::get_binary_code() const
{
	return *this;
}

void Kmer_Short::set_seq(string const &str)
{
	second_code = 0;
	for(int i = 0; i < 32; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 32; i < K - 1; i++)
	{
		second_code = (second_code << 2) + ((str[i] & 0x06) >> 1);
	}
}

string tem = "ACTG";

string Kmer_Short::get_seq() const
{
	char ch[64];
	unsigned long tmp_first_code = first_code;
	unsigned long tmp_second_code = second_code;
	int tmp_i;
	for(int i = 0; i < 32; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[31 - i] = tem[tmp_i];
	}
	for(int i = 32; i < K - 1; i++)
	{
		tmp_i = tmp_second_code & 0x0000000000000003;
		tmp_second_code = tmp_second_code >> 2;
		ch[K - i + 30] = tem[tmp_i];
	}
	return string(ch, K-1);
}

Kmer_Short &Kmer_Short::operator = (const string & str)
{
	second_code = 0;
	for(int i = 0; i < 32; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 32; i < K - 1; i++)
	{
		second_code = (second_code << 2) + ((str[i] & 0x06) >> 1);
	}
	return *this;
}

bool Kmer_Short::operator == (const Kmer_Short &tmp)
{
	if(first_code == tmp.first_code & second_code == tmp.second_code)
		return true;
	else
		return false;
}

bool Kmer_Short::operator != (const Kmer_Short &tmp)
{
	if(first_code != tmp.first_code | second_code != tmp.second_code)
		return true;
	else
		return false;
}

bool Kmer_Short::operator < (const Kmer_Short &tmp)
{
	if(first_code < tmp.first_code | (first_code == tmp.first_code & second_code < tmp.second_code))
		return true;
	else
		return false;
}

bool Kmer_Short::operator > (const Kmer_Short &tmp)
{
	if(first_code > tmp.first_code | (first_code == tmp.first_code & second_code > tmp.second_code))
		return true;
	else
		return false;
}


ostream &operator << (ostream &out, const Kmer_Short &tmp)
{
	char ch[64];
	unsigned long tmp_first_code = tmp.first_code;
	unsigned long tmp_second_code = tmp.second_code;
	int tmp_i;
	for(int i = 0; i < 32; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[31-i] = tem[tmp_i];
	}
	for(int i = 32; i < K - 1; i++)
	{
		tmp_i = tmp_second_code & 0x0000000000000003;
		tmp_second_code = tmp_second_code >> 2;
		ch[K - i + 30] = tem[tmp_i];
	}
	out << string(ch, K - 1);
	return out;
}

