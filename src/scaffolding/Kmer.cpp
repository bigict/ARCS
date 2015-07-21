#include "Kmer.h"

extern int K;

Kmer::Kmer()
{
}

Kmer::Kmer(string const &str)
{
	for(int i = 0; i < 3; i ++)
	{
		code[i] = 0;
	}
	for (int i = 0; i < K; i ++)
	{
		code[i / 32] = (code[i / 32] << 2) + ((str[i] & 0x06) >> 1);
	}

}

Kmer::Kmer(unsigned long first, unsigned long second, unsigned long third)
{
	code[0] = first;
	code[1] = second;
	code[2] = third;
}


void Kmer::set_seq(string const &str)
{
	for(int i = 0; i < 3; i ++)
	{
		code[i] = 0;
	}
	for (int i = 0; i < K; i ++)
	{
		code[i / 32] = (code[i / 32] << 2) + ((str[i] & 0x06) >> 1);
	}

}

extern string tem = "ACTG";

string Kmer::get_seq() const
{
	char ch[96];

	unsigned long tmp_code[3];
	for (int i = 0; i < 3; i++)
	{
		tmp_code[i] = code[i];
	}

	int tmp_i;
	for(int i = 0; i < 32; i++)
	{
		tmp_i = tmp_code[0] & 0x0000000000000003;
		tmp_code[0] = tmp_code[0] >> 2;
		ch[31-i] = tem[tmp_i];
	}
	for(int i = 32; i < 64; i++)
	{
		tmp_i = tmp_code[1] & 0x0000000000000003;
		tmp_code[1] = tmp_code[1] >> 2;
		ch[95 - i] = tem[tmp_i];
	}
	for(int i = 64; i < K; i++)
	{
		tmp_i = tmp_code[2] & 0x0000000000000003;
		tmp_code[2] = tmp_code[2] >> 2;
		ch[K + 63 - i] = tem[tmp_i];
	}
	return string(ch, K);

}


Kmer &Kmer::operator = (const string & str)
{
	for(int i = 0; i < 3; i ++)
	{
		code[i] = 0;
	}
	for (int i = 0; i < K; i ++)
	{
		code[i / 32] = (code[i / 32] << 2) + ((str[i] & 0x06) >> 1);
	}

}

bool Kmer::operator == (const Kmer &tmp)
{
	if(code[0] == tmp.code[0] & code[1] == tmp.code[1] & code[2] == tmp.code[2])
		return true;
	else
		return false;
}

bool Kmer::operator != (const Kmer &tmp)
{
	if(code[0] != tmp.code[0] || code[1] != tmp.code[1] || code[2] != tmp.code[2])
		return true;
	else
		return false;
}

bool Kmer::operator < (const Kmer &tmp)
{
	if(code[0] < tmp.code[0]
			|| 
			(code[0] == tmp.code[0] & code[1] < tmp.code[1])
			|| (code[0] == tmp.code[0] & code[1] == tmp.code[1] & code[2] < tmp.code[2])
	  )
		return true;
	else
		return false;
}

bool Kmer::operator > (const Kmer &tmp)
{
	if(code[0] > tmp.code[0]
			|| 
			(code[0] == tmp.code[0] & code[1] > tmp.code[1])
			|| (code[0] == tmp.code[0] & code[1] == tmp.code[1] & code[2] > tmp.code[2])
	  )
		return true;
	else
		return false;

}

/*
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
*/
