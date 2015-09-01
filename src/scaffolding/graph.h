#ifndef graph_h_
#define graph_h_

#include "component.h"

#include <iostream>
#include <list>
#include <vector>

class Graph {
public:
    struct Edge {
        size_t component_id;
        long distance;
        size_t kmer_cov; // kmer coverage
        size_t read_cov; // read coverage
        double score;
        Edge(size_t id, long d, size_t k, size_t r, double s=0.0): component_id(id), distance(d), kmer_cov(k), read_cov(r), score(s) {
        }
    };
    typedef std::list< Edge > EdgeList;
    typedef std::vector< EdgeList > NodeList;

    Graph(size_t k, size_t pair_kmer_cutoff, size_t pair_read_cutoff, double percent, size_t size, size_t genome_len) ;
    virtual ~Graph() ;

    void setPairKmerNumAndInsertSizeAndDelta(size_t pair_kmer_num, size_t insert_size, double delta) ;
    void addEdge(size_t from, size_t to, long distance, size_t kmer_num = 1, size_t read_num = 1) ;
    void scoreAndRemoveNoise(const ComponentList& components) ;
    void outputLP(std::ostream& os) ;

    size_t PAIR_KMER_NUM;
    size_t GENOME_LEN;
    size_t INSERT_SIZE;
    double DELTA;
private:
    friend std::ostream& operator<<(std::ostream& os, const Graph& g);

    NodeList _graph;
    size_t _k;
    double _percent;
    size_t _pair_read_cutoff;
    size_t _pair_kmer_cutoff;
};

#endif // graph_h_
