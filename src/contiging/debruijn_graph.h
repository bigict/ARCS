#ifndef debruijn_graph_h_
#define debruijn_graph_h_

#include "utils.h"
#include "kmer.h"

#include <iostream>
#include <string>
#include <tr1/unordered_map>

class KmerTable;

//
// https://en.wikipedia.org/wiki/De_Bruijn_graph
//
class DeBruijnGraph {
public:
    typedef std::tr1::unordered_map< Kmer, double, KmerHasher > EdgeList;

    class Node {
    public:
        Node() {
        }
        Node(const Kmer& edge, size_t n) {
            children[edge] = n;
        }
        bool empty() const {
            return children.empty() && parents.empty();
        }
        operator bool() const {
            return !empty();
        }
        size_t outdegree() const {
            return children.size();
        }
        size_t indegree() const {
            return parents.size();
        }
        
        EdgeList parents;
        EdgeList children;
    };
    typedef std::tr1::unordered_map< Kmer, Node, KmerHasher > NodeList;

    DeBruijnGraph(const KmerTable& tbl);
    virtual ~DeBruijnGraph();

    void addKmer(const Kmer& kmer, size_t weight = 1);
    void removeKmer(const Kmer& kmer);
    void removeEdge(const Kmer& node, const Kmer& edge);

    // Build a condensed graph
    void compact();
    // Remove suspicious nodes
    void removeNoise();

    const NodeList& nodelist() const {
        return _nodelist;
    }
    double average() const {
        return _average;
    }
private:
    friend std::ostream& operator << (std::ostream& os, const DeBruijnGraph& graph);

    NodeList _nodelist;
    size_t _K;
    double _average;
};

#endif // debruijn_graph_h_
