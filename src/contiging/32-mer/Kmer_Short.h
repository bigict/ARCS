#ifndef _kmer_short
#define _kmer_short

#include <string>

#include "Kmer.h"

using namespace std;

extern int K;

class Kmer;

class Kmer_Short
{
public:
	Kmer_Short();
	Kmer_Short(string const &str);
	Kmer_Short(unsigned long);
	Kmer get_next_kmer(unsigned int i) const;
	Kmer_Short get_binary_code() const;
	string get_seq() const;
	void set_seq(string const &str);
	unsigned long get_code() const;
	
	Kmer_Short & operator = (const string & str);
	//Kmer_Short & operator = (const unsigned long c_code);
	
	bool operator == (const Kmer_Short & tmp);
	bool operator < (const Kmer_Short & tmp);
	bool operator > (const Kmer_Short & tmp);
	bool operator != (const Kmer_Short & tmp);
	friend ostream& operator << (ostream& , const Kmer_Short&);

//private:
public:
	unsigned long first_code;
//	unsigned long second_code;
};

#endif
