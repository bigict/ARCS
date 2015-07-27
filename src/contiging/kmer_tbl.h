#ifndef kmer_tbl_h_
#define kmer_tbl_h_

#include <iostream>
#include <string>

#include <tr1/unordered_map>

class DeBruijnGraph;

class KmerTable {
public:
    KmerTable(size_t K);
    virtual ~KmerTable();

    bool read(std::istream& stream);
    
    void buildDeBruijn(DeBruijnGraph* graph) const;
private:
    std::tr1::unordered_map<std::string, size_t> _hash_tbl;    
    size_t _K;
};

#endif // kmer_tbl_h_
