#ifndef pair_read_h_
#define pair_read_h_

#include <string>
#include <iostream>
#include <vector>

#include "kmer_tbl.h"
#include "graph.h"

struct PairRead {
    std::string set1;
    std::string set2;
};

typedef std::vector< PairRead > PairReadList;

bool ReadPairReads(std::istream& stream1, std::istream& stream2, PairReadList& pair_reads);
bool ReadPairReads(const std::string& file1, const std::string& file2, PairReadList& pair_reads);

class InsertSizeEstimater {
public:
    InsertSizeEstimater(size_t K, size_t insert_size, const PairReadList& pair_reads, const KmerTable& kmer_tbl) : _K(K), _insert_size(insert_size), _pair_reads(pair_reads), _kmer_tbl(kmer_tbl) {
    }
    void estimate(size_t* mean, double* variance) const;
private:
    typedef std::vector< size_t > InsertSizeDistr;

    void distribution(const std::string& read1, const std::string& read2, InsertSizeDistr& insert_size_distr) const;

    size_t _K;
    size_t _insert_size;
    const PairReadList& _pair_reads;
    const KmerTable& _kmer_tbl;
};


class PairReadSet {
public:
    PairReadSet(std::istream& stream1, std::istream& stream2, size_t K, size_t INSERT_SIZE);
    PairReadSet(const std::string& file1, const std::string& file2, size_t K, size_t INSERT_SIZE);

    void buildConnectGraph(Graph& g, KmerTable& tbl, const ContigSet& contigset, const std::vector<Component>& component);
    void estimateInsertSize(const KmerTable& tbl);
    size_t size() const {
        return _pair_reads.size(); 
    }

    size_t INSERT_SIZE;
    double DELTA;
    size_t PAIR_KMER_NUM;

private:
    friend std::ostream& operator<<(std::ostream& os, const PairReadSet& p_r) ;

    void init(std::istream& stream1, std::istream& stream2);
    void estimateOnePR(const std::string& read1, const std::string& read2, const KmerTable& kmer_tbl, std::vector< int >& insert_length_list);
    std::string make_complement(std::string seq);
    void findLink(const std::string& read1, const std::string& read2, Graph& graph, const KmerTable& tbl, const ContigSet& contigset, const std::vector<Component>& components);    

    std::vector< PairRead > _pair_reads;
    size_t _k;
};

#endif // pair_read_h_
