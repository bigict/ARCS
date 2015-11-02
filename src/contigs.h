#ifndef contigs_h_
#define contigs_h_

#include <iostream>
#include <string>
#include <vector>

struct Contig {
    Contig(const std::string& seq="", size_t copy_num=0) : seq(seq), copy_num(copy_num) {
    }
    std::string seq;
    size_t copy_num;
};

typedef std::vector< Contig > ContigList;

bool ReadContigs(std::istream& stream, ContigList& contigs);
bool ReadContigs(const std::string& file, ContigList& contigs);
bool WriteContigs(std::ostream& stream, ContigList& contigs);
bool WriteContigs(const std::string& file, ContigList& contigs);
size_t GenomeLength(const ContigList& contigs, size_t K);

// FA(contig file) Format Specification: file_name(.fa)
// Syntax
//    <seqname, copy_num>    :=    [A-Za-z0-9_.:-]+  [0-9]+
//    <seq>    :=    [A-Za-z\n\.~]+
// Requirements
//    1. The <seqname> appears right after '>', and seqname and copy_num in the small line. 

class ContigReader {
public:
    ContigReader(std::istream& stream) : _stream(stream) {
    }
    bool read(Contig& contig);

private:
    std::istream& _stream;
};

#endif // contigs_h_
