#ifndef pair_read_h_
#define pair_read_h_

#include <iostream>
#include <string>
#include <vector>

#include "component.h"
#include "kmer_tbl.h"

class KmerTable;
class Graph;

struct PairRead {
    std::string set1;
    std::string set2;
};

typedef std::vector< PairRead > PairReadList;
bool ReadPairReads(std::istream& stream1, std::istream& stream2, PairReadList& pair_reads);
bool ReadPairReads(const std::string& file1, const std::string& file2, PairReadList& pair_reads);

class InsertSizeEstimater {
public:
    InsertSizeEstimater(size_t K, size_t insert_size, const PairReadList& pair_reads, const KmerList& hash_tbl) : _K(K), _insert_size(insert_size), _pair_reads(pair_reads), _hash_tbl(hash_tbl) {
    }

    void estimate(size_t* insert_size, double* delta);
private:
    typedef std::vector< int > InsertSizeDistr;

    void estimateOnePR(const std::string& read1, const std::string& read2, InsertSizeDistr& insert_size_distr);

    size_t _K;
    size_t _insert_size;
    const PairReadList& _pair_reads;
    const KmerList& _hash_tbl;
};

class ConnectGraphBuilder {
public:
    ConnectGraphBuilder(size_t K, size_t insert_size, const PairReadList& pair_reads, const KmerList& hash_tbl, const ComponentList& components) : _K(K), _insert_size(insert_size), _pair_reads(pair_reads), _hash_tbl(hash_tbl), _components(components) {
    }

    size_t build(Graph* g) const;
private:
    size_t addEdge(const std::string& read1, const std::string& read2, Graph* graph) const; 

    size_t _K;
    size_t _insert_size;
    const PairReadList& _pair_reads;
    const KmerList& _hash_tbl;
    const ComponentList& _components;
};

#endif // pair_read_h_
