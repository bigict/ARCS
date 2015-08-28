#ifndef contigs_h_
#define contigs_h_

#include <string>
#include <iostream>
#include <vector>

struct Contig {
    std::string contig;
    size_t copy_num;
    Contig(const std::string contig="", size_t copy_num=0) : contig(contig), copy_num(copy_num) {
    }
};

class ContigSet {
public:
    ContigSet(std::istream& stream, size_t K); 
    bool read(Contig& contig);
    friend std::ostream& operator<<(std::ostream& os, const ContigSet& contigs);

    std::vector< Contig > _contigs;
    size_t GENOME_LEN;

private:
    std::istream& _stream;
    size_t _k;
};


// FA(contig file) Format Specification: file_name(.fa)
// Syntax
//    <seqname, copy_num>    :=    [A-Za-z0-9_.:-]+  [0-9]+
//    <seq>    :=    [A-Za-z\n\.~]+
// Requirements
//    1. The <seqname> appears right after '>', and seqname and copy_num in the small line. 


#endif // contigs_h_
