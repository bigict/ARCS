#ifndef condensed_debruijn_graph_h_
#define condensed_debruijn_graph_h_

#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

//
// https://en.wikipedia.org/wiki/De_Bruijn_graph
//

class CondensedDeBruijnGraph {
public:
    struct CondensedEdge {
        CondensedEdge() {
        }
        CondensedEdge(const std::string& seq, size_t copy_num) : seq(seq), copy_num(copy_num) {
        }

        std::string seq;
        size_t copy_num;

        size_t length() const {
            return seq.length();
        }
    };

    CondensedDeBruijnGraph(size_t K) : _K(K) {
    }

    void addEdge(const std::string& seq, size_t copy_num);
    void removeEdge(const std::string& seq);

    typedef std::vector< CondensedEdge > CondensedEdgeList;
    CondensedEdgeList _indexer;
    typedef std::map< char, size_t > EdgeList;
    typedef std::unordered_map< std::string, EdgeList >  NodeList;
    NodeList _parents;
    NodeList _children;
private:
    friend std::ostream& operator<<(std::ostream& os, const CondensedDeBruijnGraph& g);

    size_t _K;
};

#endif // condensed_debruijn_graph_h_
