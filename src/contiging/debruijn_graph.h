#ifndef debruijn_graph_h_
#define debruijn_graph_h_

#include "utils.h"
#include "kmer.h"

#include <string>
#include <tr1/unordered_map>

class KmerTable;

//
// https://en.wikipedia.org/wiki/De_Bruijn_graph
//
class DeBruijnGraph {
public:
    DeBruijnGraph(const KmerTable& tbl);
    virtual ~DeBruijnGraph();

    void addKmer(const Kmer& kmer, size_t weight = 1);
    void removeKmer(const Kmer& kmer);

    // Build a condensed graph
    void compact();
    
    // Remove suspicious nodes
    void removeNoise();
private:
    typedef std::tr1::unordered_map< Kmer, size_t, KmerHasher > EdgeList;

    class Node {
    public:
        Node() {
        }
        Node(const Kmer& edge, size_t n) {
            children[edge] = n;
        }
        bool empty() const {
            return children.empty();
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
        size_t average() const;
        
        EdgeList parents;
        EdgeList children;
    };

    struct NodeKey {
        NodeKey(size_t K) : _K(K) {}
        Kmer operator()(const Kmer& kmer) const {
            size_t k = std::min(kmer.length(), _K);
            return kmer.subKmer(kmer.length() -  k);
        }
    private:
        size_t _K;
    };

    typedef std::tr1::unordered_map< Kmer, Node, KmerHasher > NodeList;
    NodeList _nodelist;
    size_t _K;
    double _average;
};

#endif // debruijn_graph_h_
