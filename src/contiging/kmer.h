#ifndef kmer_h_
#define kmer_h_

#include <string>

#include <boost/multiprecision/cpp_int.hpp>

class Nucleotide {
public:
    enum {
        Adenine  = 0x00, 
        Cytosine = 0x01, 
        Guanine  = 0x10, 
        Thymine  = 0x11
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

class Kmer {
public:
    Kmer();
    Kmer(const std::string& sequence);
    Kmer(const std::string& sequence, size_t i, size_t j);
    Kmer(const Kmer& o);
    virtual ~Kmer();

    Kmer subKmer(size_t i, size_t j) const;

    const std::string sequence() const;
    void sequence(const std::string& seq) {
        sequence(seq, 0, std::string::npos);
    }
    void sequence(const std::string& seq, size_t i, size_t j);
    size_t length() const { return _length; }

    Kmer& operator = (const std::string& sequence);
    Kmer& operator = (const Kmer& o);
    bool operator == (const Kmer& o) const;
    bool operator < (const Kmer& o) const;
    bool operator > (const Kmer& o) const;
    bool operator != (const Kmer& o) const;

private:
    boost::multiprecision::cpp_int _data;
    size_t _length;
};

#endif // kmer_h_
