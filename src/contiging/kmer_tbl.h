#ifndef kmer_tbl_h_
#define kmer_tbl_h_

#include "kmer.h"

#include <iostream>
#include <string>
#include <tr1/unordered_map>
#include <vector>

class DNASeq;
class DeBruijnGraph;

class KmerTable {
public:
    KmerTable(size_t K, double avg_quality=0, double min_quality=0, double percent=1.0, size_t read_cutoff=-1, bool do_reversed=true);
    virtual ~KmerTable();

    size_t K() const { return _K; }

    bool read(std::istream& stream);
    bool read(const std::string& stream);
    bool read(const std::vector< std::string >& stream_list);
    bool write(std::ostream& stream) const;
    void buildDeBruijn(DeBruijnGraph* graph) const;

    void statistics(double* average, double* variance) const;
private:
    void addRead(const DNASeq& read);
    bool isValid(const DNASeq& read, size_t i, size_t j) const;

    typedef std::tr1::unordered_map< Kmer, size_t, KmerHasher > KmerList;
    KmerList _hash_tbl;    
    size_t _K;
    double _avg_quality;
    double _min_quality;
    double _percent;
    size_t _read_cutoff;
    bool _do_reversed;
};

#endif // kmer_tbl_h_
