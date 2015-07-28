#ifndef kmer_tbl_h_
#define kmer_tbl_h_

#include "kmer.h"

#include <iostream>
#include <map>
#include <string>
#include <tr1/unordered_map>

class DeBruijnGraph;

class KmerTable {
public:
    KmerTable(size_t K, bool do_filter=false);
    virtual ~KmerTable();

    size_t K() const { return _K; }

    bool read(std::istream& stream);
    void buildDeBruijn(DeBruijnGraph* graph) const;

private:
    typedef std::tr1::unordered_map< Kmer, size_t, KmerHasher > KmerList;
    KmerList _hash_tbl;    
    size_t _K;
    size_t _do_filter;
};

#endif // kmer_tbl_h_
