#include "Kmer.h"

extern int K;
Kmer::Kmer()
{
}

unsigned long Kmer::get_code()const{
	return first_code;
}

Kmer::Kmer(string const &str)
{
	first_code = 0;
	for(int i = 0; i < K; i++)
	{
		if(str[i] == 'N')
			first_code = (first_code << 2) | (rand() & 0x3);
		else
			first_code = (first_code << 2) | ((str[i] & 0x06) >> 1);
	}
}

Kmer::Kmer(unsigned long first)
{
	first_code = first;
	// cout << first << " " << first_code << endl;
}

Kmer Kmer::get_binary_code() const
{
	return *this;
}

void Kmer::set_seq(string const &str)
{
	first_code = 0;
	for(int i = 0; i < K; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
}

extern string tem;

string Kmer::get_seq() const
{
	char ch[32];
	unsigned long tmp_first_code = first_code;
	int tmp_i;
	for(int i = 0; i < K; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[K - i -1] = tem[tmp_i];
	}
	return string(ch, K);

}

Kmer_Short Kmer::get_next_node()
{
	unsigned long tmp_first_code = first_code;
	tmp_first_code = tmp_first_code << (64 - 2*K + 2);
	tmp_first_code = tmp_first_code >> (64 - 2*K + 2);
	return Kmer_Short(tmp_first_code);
}

Kmer_Short Kmer::get_pre_node()
{
	unsigned long tmp_first_code = first_code;
	tmp_first_code = tmp_first_code >> 2;
	return Kmer_Short(tmp_first_code);
}

Kmer &Kmer::operator = (const string & str)
{
	first_code = 0;
	for(int i = 0; i < K; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	return *this;
}

bool Kmer::operator == (const Kmer &tmp)
{
	if(first_code == tmp.first_code)
		return true;
	else
		return false;
}

bool Kmer::operator != (const Kmer &tmp)
{
	if(first_code != tmp.first_code)
		return true;
	else
		return false;
}

bool Kmer::operator < (const Kmer &tmp)
{
	if(first_code < tmp.first_code )
		return true;
	else
		return false;
}

bool Kmer::operator > (const Kmer &tmp)
{
	if(first_code > tmp.first_code)
		return true;
	else
		return false;
}

ostream &operator << (ostream &out, const Kmer &tmp)
{
	char ch[32];
	unsigned long tmp_first_code = tmp.first_code;
	int tmp_i;
	for(int i = 0; i < K; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[K-i-1] = tem[tmp_i];
	}
	out << string(ch, K);
	return out;
}
