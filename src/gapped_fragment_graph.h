#ifndef gapped_fragment_graph_h_
#define gapped_fragment_graph_h_

#include "component.h"

#include <iostream>
#include <list>
#include <vector>
#include <set>
#include <map>

class GappedFragmentGraph {
public:
    struct Edge {
        size_t component_id;
        long distance;
        std::vector<double> distances;
        size_t kmer_cov; // kmer coverage
        size_t read_cov; // read coverage
        double score;
        Edge(size_t id, long d, size_t k, size_t r, double s=0.0): component_id(id), distance(d), kmer_cov(k), read_cov(r), score(s) {
            distances.push_back( static_cast<double>(d) );
        }
    };
    typedef std::list< Edge > EdgeList;
    typedef std::vector< EdgeList > NodeList;
    typedef std::set< size_t > RepeateList;
    typedef std::set< std::pair<size_t, size_t> > ReverseMapEle;
    typedef std::map<size_t, ReverseMapEle> ReverseMap; // node to <node, pair_read_num>

    GappedFragmentGraph(size_t K, size_t pair_kmer_cutoff, size_t pair_read_cutoff, double percent, size_t size, size_t genome_len) ;
    virtual ~GappedFragmentGraph() ;

    void addEdge(size_t from, size_t to, long distance, size_t kmer_num = 1, size_t read_num = 1) ;
    void scoreAndRemoveNoise(const ComponentList& components, size_t read_len) ;
    bool deepMoreThan2(size_t startComp) ;
    bool reverseDeepMoreThan2(size_t startComp, ReverseMap &reverseMap) ;
    void outputLP(std::ostream& os, const ComponentList& components) ;

    size_t PAIR_KMER_NUM;
    size_t GENOME_LEN;
    size_t INSERT_SIZE;
    double DELTA;
    std::map<size_t, int> node2PairReads;
private:
    friend std::ostream& operator<<(std::ostream& os, const GappedFragmentGraph& g);

    NodeList _nodelist;
    RepeateList _repeateList;
    size_t _K;
    double _percent;
    size_t _pair_read_cutoff;
    size_t _pair_kmer_cutoff;
};

#endif // gapped_fragment_graph_h_
