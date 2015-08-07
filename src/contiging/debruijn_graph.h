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
    typedef std::tr1::unordered_map< Kmer, size_t, KmerHasher > EdgeList;

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

    void outputCdbgGraphAndComponent0(const std::string& cdbg_file, const std::string& component0_file) const;

    struct NodeKey {
        NodeKey(size_t K) : _K(K) {}
        Kmer operator()(const Kmer& kmer) const {
            size_t k = std::min(kmer.length(), _K);
            return kmer.subKmer(kmer.length() -  k);
        }
    private:
        size_t _K;
    };

private:
    friend std::ostream& operator << (std::ostream& os, const DeBruijnGraph& graph);

    NodeList _nodelist;
    size_t _K;
    double _average;
};

#endif // debruijn_graph_h_
