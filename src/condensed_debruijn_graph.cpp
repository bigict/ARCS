#include "condensed_debruijn_graph.h"

#include <boost/assert.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.CondensedDeBruijnGraph"));

void CondensedDeBruijnGraph::addEdge(const std::string& seq, size_t copy_num) {
    BOOST_ASSERT(seq.length() >= _K);

    CondensedEdge edge(seq, copy_num);
    _indexer.push_back(edge);
    size_t idx = _indexer.size() - 1;

    // parents
    {
        std::string prefix = edge.seq.substr(0, _K - 1);
        char next = edge.seq[_K - 1];
        NodeList::iterator  i = _parents.find(prefix);
        if (i != _parents.end()) {
            i->second[next] = idx;
        } else {
            EdgeList parent = boost::assign::map_list_of(next, idx);
            _parents[prefix] = parent;
        }
    }
    // children
    {
        std::string suffix = edge.seq.substr(edge.seq.length() - (_K - 1), _K - 1);
        char prev = edge.seq[edge.seq.length() - _K];
        NodeList::iterator  i = _children.find(suffix);
        if (i != _children.end()) {
            i->second[prev] = idx;
        } else {
            EdgeList child = boost::assign::map_list_of(prev, idx);
            _children[suffix] = child;
        }
    }
}

std::ostream& operator<<(std::ostream& os, const CondensedDeBruijnGraph& g) {
    os << boost::format("edges=[%d], parents=[%d], children=[%d]") % g._indexer.size() % g._parents.size() % g._children.size();
    return os;
}
