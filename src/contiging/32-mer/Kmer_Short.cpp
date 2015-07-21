#include "Kmer_Short.h"

Kmer_Short::Kmer_Short()
{
}

Kmer_Short::Kmer_Short(string const &str)
{
	first_code = 0;
	for(int i = 0; i < K - 1; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
}

Kmer_Short::Kmer_Short(unsigned long first)
{
	first_code = first;
}

unsigned long Kmer_Short::get_code()const{
	return first_code;
}

Kmer Kmer_Short::get_next_kmer(unsigned int i) const
{
	unsigned long first = first_code;
	first = first << 2;
	first += i;
	return Kmer(first);
}

Kmer_Short Kmer_Short::get_binary_code() const
{
	return *this;
}

void Kmer_Short::set_seq(string const &str)
{
	first_code = 0;
	for(int i = 0; i < K - 1; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
}

string tem = "ACTG";

string Kmer_Short::get_seq() const
{
	char ch[32];
	unsigned long tmp_first_code = first_code;
	int tmp_i;
	for(int i = 0; i < K - 1; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[K-i-2] = tem[tmp_i];
	}
	return string(ch, K-1);
}

Kmer_Short &Kmer_Short::operator = (const string & str)
{
	first_code = 0;
	for(int i = 0; i < K - 1; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	return *this;
}

bool Kmer_Short::operator == (const Kmer_Short &tmp)
{
	if(first_code == tmp.first_code)
		return true;
	else
		return false;
}

bool Kmer_Short::operator != (const Kmer_Short &tmp)
{
	if(first_code != tmp.first_code)
		return true;
	else
		return false;
}

bool Kmer_Short::operator < (const Kmer_Short &tmp)
{
	if(first_code < tmp.first_code)
		return true;
	else
		return false;
}

bool Kmer_Short::operator > (const Kmer_Short &tmp)
{
	if(first_code > tmp.first_code)
		return true;
	else
		return false;
}


ostream &operator << (ostream &out, const Kmer_Short &tmp)
{
	char ch[32];
	unsigned long tmp_first_code = tmp.first_code;
	int tmp_i;
	for(int i = 0; i < K - 1; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[K - 2 - i] = tem[tmp_i];
	}
	out << string(ch, K - 1);
	return out;
}

