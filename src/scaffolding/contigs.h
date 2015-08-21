#ifndef contigs_h_
#define contigs_h_

#include <iostream>
#include <string>
#include <vector>

struct Contig {
    Contig(const std::string contig="", size_t copy_num=0) : contig(contig), copy_num(copy_num) {
    }
    size_t id;
    std::string contig;
    size_t copy_num;
};

class ContigSet {
public:
    typedef std::vector< Contig > ContigList;

    ContigSet(std::istream& stream, size_t K); 
    ContigSet(const std::string& filename, size_t K); 

    ContigList contigs;
    size_t GENOME_LEN;

private:
    void init(std::istream& stream);

    friend std::ostream& operator<<(std::ostream& os, const ContigSet& contigs);
    size_t _K;
};

class ContigReader {
public:
    ContigReader(std::istream& stream) : _stream(stream) {
    }
    bool read(Contig& contig);

private:
    std::istream& _stream;
};

#endif // contigs_h_
