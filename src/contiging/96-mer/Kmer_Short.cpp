#include "Kmer_Short.h"

Kmer_Short::Kmer_Short()
{
}

Kmer_Short::Kmer_Short(string const &str)
{
	third_code = 0;
	for(int i = 0; i < 32; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 32; i < 64; i++)
	{
		second_code = (second_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 64; i < K - 1; i++)
	{
		third_code = (third_code << 2) + ((str[i] & 0x06) >> 1);
	}
}

Kmer_Short::Kmer_Short(unsigned long first, unsigned long second, unsigned long third)
{
	first_code = first;
	second_code = second;
	third_code = third;
}

unsigned long Kmer_Short::get_code() const{
	return first_code | second_code | third_code;
}

Kmer Kmer_Short::get_next_kmer(unsigned int i) const
{
	unsigned long third = third_code;

	third = third << 2;
	third += i;
	return Kmer(first_code, second_code, third);
}

Kmer_Short Kmer_Short::get_binary_code() const
{
	return *this;
}

Kmer_Short Kmer_Short::get_next_short(unsigned int i) const
{
	unsigned long third = third_code << 2;
	third |= i;
	return Kmer(first_code, second_code, third).get_next_node();
}

Kmer_Short Kmer_Short::get_next_short() const
{
	return Kmer(first_code, second_code, third_code).get_next_node();
}

void Kmer_Short::set_seq(string const &str)
{
	third_code = 0;
	for(int i = 0; i < 32; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 32; i < 64; i++)
	{
		second_code = (second_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 64; i < K - 1; i++)
	{
		third_code = (third_code << 2) + ((str[i] & 0x06) >> 1);
	}

}

string tem = "ACTG";

string Kmer_Short::get_seq() const
{
	char ch[96];
	unsigned long tmp_first_code = first_code;
	unsigned long tmp_second_code = second_code;
	unsigned long tmp_third_code = third_code;
	int tmp_i;
	for(int i = 0; i < 32; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[31 - i] = tem[tmp_i];
	}
	for(int i = 32; i < 64; i++)
	{
		tmp_i = tmp_second_code & 0x0000000000000003;
		tmp_second_code = tmp_second_code >> 2;
		ch[95 - i] = tem[tmp_i];
	}
	for(int i = 64; i < K - 1; i++)
	{
		tmp_i = tmp_third_code & 0x0000000000000003;
		tmp_third_code = tmp_third_code >> 2;
		ch[K + 62 - i] = tem[tmp_i];
	}
	return string(ch, K-1);
}

Kmer_Short &Kmer_Short::operator = (const string & str)
{
	third_code = 0;
	for(int i = 0; i < 32; i++)
	{
		first_code = (first_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 32; i < 64; i++)
	{
		second_code = (second_code << 2) + ((str[i] & 0x06) >> 1);
	}
	for(int i = 64; i < K - 1; i++)
	{
		third_code = (third_code << 2) + ((str[i] & 0x06) >> 1);
	}
	return *this;
}
bool Kmer_Short::operator == (const Kmer_Short &tmp)
{
	if(first_code == tmp.first_code & second_code == tmp.second_code & third_code == tmp.third_code)
		return true;
	else
		return false;
}

bool Kmer_Short::operator != (const Kmer_Short &tmp)
{
	if(first_code != tmp.first_code || second_code != tmp.second_code || third_code != tmp.third_code)
		return true;
	else
		return false;
}

bool Kmer_Short::operator < (const Kmer_Short &tmp)
{
	if(first_code < tmp.first_code 
			|| 
			(first_code == tmp.first_code & second_code < tmp.second_code)
			||
			(first_code == tmp.first_code & second_code == tmp.second_code & third_code < tmp.third_code)
	  )
		return true;
	else
		return false;
}

bool Kmer_Short::operator > (const Kmer_Short &tmp)
{
	if(first_code > tmp.first_code 
			|| 
			(first_code == tmp.first_code & second_code > tmp.second_code)
			||
			(first_code == tmp.first_code & second_code == tmp.second_code & third_code > tmp.third_code)
	  )
		return true;
	else
		return false;
}


ostream &operator << (ostream &out, const Kmer_Short &tmp)
{
	char ch[96];
	unsigned long tmp_first_code = tmp.first_code;
	unsigned long tmp_second_code = tmp.second_code;
	unsigned long tmp_third_code = tmp.third_code;
	int tmp_i;
	for(int i = 0; i < 32; i++)
	{
		tmp_i = tmp_first_code & 0x0000000000000003;
		tmp_first_code = tmp_first_code >> 2;
		ch[31-i] = tem[tmp_i];
	}
	for(int i = 32; i < 64; i++)
	{
		tmp_i = tmp_second_code & 0x0000000000000003;
		tmp_second_code = tmp_second_code >> 2;
		ch[95 - i] = tem[tmp_i];
	}
	for(int i = 64; i < K-1; i++)
	{
		tmp_i = tmp_third_code & 0x0000000000000003;
		tmp_third_code = tmp_third_code >> 2;
		ch[K + 62 - i] = tem[tmp_i];
	}

	out << string(ch, K - 1);
	return out;
}

