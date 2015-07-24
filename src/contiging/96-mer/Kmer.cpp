#include "Kmer.h"

extern int K;

Kmer::Kmer()
{
}

Kmer::Kmer(string const &str)
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
	for(int i = 64; i < K; i ++)
	{
		third_code = (third_code << 2) + ((str[i] & 0x06) >> 1);
	}

}

Kmer::Kmer(unsigned long first, unsigned long second, unsigned long third)
{
	first_code = first;
	second_code = second;
	third_code = third;
}

unsigned long Kmer::get_code() const{
	return first_code | second_code | third_code;
}

Kmer Kmer::get_binary_code() const
{
	return *this;
}

void Kmer::set_seq(string const &str)
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
	for(int i = 64; i < K; i ++)
	{
		third_code = (third_code << 2) + ((str[i] & 0x06) >> 1);
	}

}

extern string tem;

string Kmer::get_seq() const
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
		ch[31-i] = tem[tmp_i];
	}
	for(int i = 32; i < 64; i++)
	{
		tmp_i = tmp_second_code & 0x0000000000000003;
		tmp_second_code = tmp_second_code >> 2;
		ch[95 - i] = tem[tmp_i];
	}
	for(int i = 64; i < K; i++)
	{
		tmp_i = tmp_third_code & 0x0000000000000003;
		tmp_third_code = tmp_third_code >> 2;
		ch[K + 63 - i] = tem[tmp_i];
	}
	return string(ch, K);

}

Kmer_Short Kmer::get_next_node()
{
	unsigned long tmp_first_code = first_code;
	unsigned long tmp_second_code = second_code;
	unsigned long tmp_third_code = third_code;
	tmp_first_code = tmp_first_code << 2;
	unsigned long ff = 0xc000000000000000;
	unsigned long last = ff & tmp_second_code;
	last = last >> 62;
	tmp_first_code += last;
	tmp_second_code = tmp_second_code << 2;
	ff = 0x0000000000000003;
	ff = ff << (2*(K - 65));
	last = ff & tmp_third_code;
	last = last >> (2*(K - 65));
	ff = ~ff;
	tmp_third_code &= ff;
	tmp_second_code += last;
	return Kmer_Short(tmp_first_code, tmp_second_code, tmp_third_code);
}

Kmer_Short Kmer::get_pre_node()
{
	unsigned long tmp_third_code = third_code;
	tmp_third_code = tmp_third_code >> 2;
	return Kmer_Short(first_code, second_code, tmp_third_code);
}

Kmer &Kmer::operator = (const string & str)
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
	for(int i = 64; i < K; i++)
	{
		third_code = (third_code << 2) + ((str[i] & 0x06) >> 1);
	}
	return *this;
}

bool Kmer::operator == (const Kmer &tmp)
{
	if(first_code == tmp.first_code & second_code == tmp.second_code & third_code == tmp.third_code)
		return true;
	else
		return false;
}

bool Kmer::operator != (const Kmer &tmp)
{
	if(first_code != tmp.first_code || second_code != tmp.second_code || third_code != tmp.third_code)
		return true;
	else
		return false;
}

bool Kmer::operator < (const Kmer &tmp)
{
	if(first_code < tmp.first_code
			|| 
			(first_code == tmp.first_code & second_code < tmp.second_code)
			|| (first_code == tmp.first_code & second_code == tmp.second_code & third_code < tmp.third_code)
	  )
		return true;
	else
		return false;
}

bool Kmer::operator > (const Kmer &tmp)
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

ostream &operator << (ostream &out, const Kmer &tmp)
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
	for(int i = 64; i < K; i++)
	{
		tmp_i = tmp_third_code & 0x0000000000000003;
		tmp_third_code = tmp_third_code >> 2;
		ch[K + 63 - i] = tem[tmp_i];
	}
	out << string(ch, K);
	return out;
}
