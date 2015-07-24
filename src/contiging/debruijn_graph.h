#ifndef debruijn_graph_h_
#define debruijn_graph_h_

#include <string>
#include <tr1/unordered_map>

class DeBruijnGraph {
public:
    DeBruijnGraph();
    virtual ~DeBruijnGraph();

private:
    class Node {
    public:
        size_t count[4];
    };

    std::tr1::unordered_map<std::string, size_t> _hash_tbl;    
};

#endif // debruijn_graph_h_
