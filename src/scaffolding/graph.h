#ifndef graph_h_
#define graph_h_

#include "kmer.h"
#include "component.h"

#include <iostream>
#include <vector>
#include <list>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
struct Edge {
    size_t component_id;
    long dis;
    size_t kmer_cov;
    size_t read_cov;
    double score;
    Edge(size_t id, long d, size_t k, size_t r, double s=0.0): component_id(id), dis(d), kmer_cov(k), read_cov(r), score(s) {
	}
};

class Graph {
public:
    Graph(size_t k, size_t pair_kmer_cutoff, size_t pair_read_cutoff, double percent, size_t size, size_t genome_len) ;
    virtual ~Graph() ;

    void setPairKmerNumAndInsertSizeAndDelta(size_t pair_kmer_num, size_t insert_size, double delta) ;
    void addEdge(size_t from, size_t to, long dis, bool isread) ;
    void scoreAndRemoveNoise(const std::vector<Component>& com) ;
    void outputLP(std::ostream& os) ;
    friend std::ostream& operator<<(std::ostream& os, const Graph& g) ;

private:
    double score(size_t leni, size_t lenj, long d, size_t c) ;
    void initialize_gaussian() ;
    double get_Gaussian(int i) ;
    double log_poisson(double lam, int i) ;
    //double get_Gaussian(int i) ;
public:
    typedef std::list<Edge> GraphEdge;
    typedef std::vector< GraphEdge > GraphNode;
private:
    GraphNode _graph;
    size_t _k;
    double _percent;
    size_t _pair_read_cutoff;
    size_t _pair_kmer_cutoff;
    size_t PAIR_KMER_NUM;
    size_t GENOME_LEN;
    size_t INSERT_SIZE;
    double DELTA;
    std::vector<double> Gaussian;
};

#endif
