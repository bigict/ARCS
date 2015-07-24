#ifndef kmer_h_
#define kmer_h_

#include <string>

#include <boost/multiprecision/cpp_int.hpp>

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

    Kmer& operator = (const std::string& sequence);
    bool operator == (const Kmer& o);
    bool operator < (const Kmer& o);
    bool operator > (const Kmer& o);
    bool operator != (const Kmer& o);

private:
    boost::multiprecision::cpp_int _data;
};

#endif // kmer_h_
