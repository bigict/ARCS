#ifndef _kmer
#define _kmer

#include <string>
#include <iostream>

using namespace std;

extern int K;

class Kmer
{
public:
	Kmer();
	Kmer(string const &str);
	Kmer(unsigned long, unsigned long, unsigned long);
	string get_seq() const;
	void set_seq(string const &str);
	
	Kmer & operator = (const string & str);
	bool operator == (const Kmer & tmp);
	bool operator < (const Kmer & tmp);
	bool operator > (const Kmer & tmp);
	bool operator != (const Kmer & tmp);
	
//	friend ostream& operator<< (ostream& , const Kmer&);
	
public:
	unsigned long code[3];
};

#endif
