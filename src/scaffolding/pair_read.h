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

class PairReadSet {
public:
    PairReadSet(std::istream& stream1, std::istream& stream2, size_t K, size_t INSERT_SIZE);
    PairReadSet(const std::string& file1, const std::string& file2, size_t K, size_t INSERT_SIZE);

    void buildConnectGraph(Graph& g, KmerTable& tbl, const std::vector<Component>& component);
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
    void findLink(const std::string& read1, const std::string& read2, Graph& graph, const KmerTable& tbl, const std::vector<Component>& components);    

    std::vector< PairRead > _pair_reads;
    size_t _k;
};

#endif // pair_read_h_
