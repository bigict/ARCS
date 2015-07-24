#ifndef kmer_tbl_h_
#define kmer_tbl_h_

#include <string>

#include <tr1/unordered_map>

class KmerTable {
public:
    KmerTable();
    virtual ~KmerTable();

private:
    std::tr1::unordered_map<std::string, size_t> _hash_tbl;    
};

#endif // kmer_tbl_h_
