#ifndef _kmer
#define _kmer

#include <string>
#include <iostream>

#include "Kmer_Short.h"

using namespace std;

extern int K;

class Kmer_Short;

class Kmer
{
public:
	Kmer();
	Kmer(string const &str);
	Kmer(unsigned long, unsigned long);
	Kmer get_binary_code() const;
	string get_seq() const;
	unsigned long get_code()const;
	void set_seq(string const &str);
	
	Kmer & operator = (const string & str);
//	Kmer & operator = (const unsigned long c_code);
	bool operator == (const Kmer & tmp);
	bool operator < (const Kmer & tmp);
	bool operator > (const Kmer & tmp);
	bool operator != (const Kmer & tmp);
	
	friend ostream& operator<< (ostream& , const Kmer&);
	
	Kmer_Short get_next_node();
	Kmer_Short get_pre_node();
//private:
public:
	unsigned long first_code;
	unsigned long second_code;
};

#endif
