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
private:
    class Node {
    public:
        Node() {
            memset(count, 0, _countof(count));
            memset(length, 0, _countof(length));
        }
        Node(Nucleotide::Code nucleotide, size_t n) {
            memset(count, 0, _countof(count));
            memset(length, 0, _countof(length));
            count[nucleotide] = n;
        }
        Node(char nucleotide, size_t n) {
            memset(count, 0, _countof(count));
            memset(length, 0, _countof(length));
            count[Nucleotide::char2code(nucleotide)] = n;
        }
        bool empty() const {
            for (size_t i = 0; i < _countof(count); ++i) {
                if (count[i] > 0) {
                    return false;
                }
            }
            return true;
        }
        operator bool() const {
            return !empty();
        }
        size_t outdegree() const {
            size_t n = 0;
            for (size_t i = 0; i < _countof(count); ++i) {
                if (count[i] > 0) {
                    ++n;
                }
            }
            return n;
        }

        size_t count[4];
        size_t length[4];
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
};

#endif // debruijn_graph_h_
