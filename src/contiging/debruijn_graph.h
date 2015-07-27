#ifndef debruijn_graph_h_
#define debruijn_graph_h_

#include "kmer.h"

#include <string>
#include <tr1/unordered_map>

class KmerTable;

class DeBruijnGraph {
public:
    DeBruijnGraph(const KmerTable& tbl);
    virtual ~DeBruijnGraph();

    void addKmer(const Kmer& kmer, size_t weight = 1);
    void removeKmer(const Kmer& kmer);
private:
    class Node {
    public:
        Node() {
            memset(count, 0, sizeof(count) / sizeof(count[0]));
        }
        Node(Nucleotide::Code code, size_t n) {
            memset(count, 0, sizeof(count) / sizeof(count[0]));
            count[code] = n;
        }
        Node(char nucleotide, size_t n) {
            memset(count, 0, sizeof(count) / sizeof(count[0]));
            count[Nucleotide::char2code(nucleotide)] = n;
        }
        bool empty() const {
            for (size_t i = 0; i < sizeof(count) / sizeof(count[0]); ++i) {
                if (count[i] > 0) {
                    return false;
                }
            }
            return true;
        }
        operator bool() const {
            return !empty();
        }
        size_t count[4];
    };

    typedef std::tr1::unordered_map< Kmer, Node, KmerHasher > NodeList;
    //std::tr1::unordered_map< Kmer, Node, KmerHasher> _hash_tbl;    
    NodeList _nodelist;
};

#endif // debruijn_graph_h_
