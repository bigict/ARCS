#ifndef kmer_h_
#define kmer_h_

#include <string>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>

typedef boost::multiprecision::gmp_int bigint;

//
// Utilities for encodeing Nucleotide
//
class Nucleotide {
public:
    enum Code {
        Adenine  = 0x00, 
        Cytosine = 0x01, 
        Guanine  = 0x02, 
        Thymine  = 0x03,
        NUM      = 0x04
    };
    
    static int char2code(char chr) {
        chr = std::toupper(chr);
        switch (chr) {
        case 'A': return Adenine;
        case 'C': return Cytosine;
        case 'G': return Guanine;
        case 'T': return Thymine;
        }
        return Adenine;
    }
    
    static char code2char(int code) {
        switch (code) {
        case Adenine:  return 'A';
        case Cytosine: return 'C';
        case Guanine:  return 'G';
        case Thymine:  return 'T';
        }
        return 'A';
    }
};

//
// A node in DeBruij graph. Here we encode the nucleotide to reduce
// memeory. You can restore sequence by calling the sequence member 
// function without losing any information.
//
class Kmer {
public:
    Kmer();
    Kmer(const std::string& sequence);
    Kmer(const std::string& sequence, size_t i, size_t j=std::string::npos);
    Kmer(const Kmer& o);
    virtual ~Kmer();

    Kmer subKmer(size_t i, size_t j=-1) const;

    Nucleotide::Code nucleotide(size_t i) const;
    const std::string sequence() const;
    void sequence(const std::string& seq) {
        sequence(seq, 0, std::string::npos);
    }
    void sequence(const std::string& seq, size_t i, size_t j=std::string::npos);
    size_t length() const { return _length; }

    Kmer& operator = (const std::string& sequence);
    Kmer& operator = (const Kmer& o);
    Kmer operator + (const char c);
    Kmer operator + (const std::string& sequence);
    Kmer operator + (const Kmer& o);
    Kmer& operator += (const char c);
    Kmer& operator += (const std::string& sequence);
    Kmer& operator += (const Kmer& o);
    bool operator < (const Kmer& o) const { return _data < o._data; }
    bool operator > (const Kmer& o) const { return _data > o._data; }
    bool operator == (const Kmer& o) const;
    bool operator != (const Kmer& o) const;

    size_t hash() const;
private:
    friend std::ostream& operator << (std::ostream& os, const Kmer& kmer);

    bigint _data;
    size_t _length;
};

struct KmerHasher {
        std::size_t operator()(const Kmer& o) const {
            return o.hash();
        }
};

#endif // kmer_h_
