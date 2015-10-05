#ifndef debruijn_graph_h_
#define debruijn_graph_h_

#include "arff.h"
#include "kmer.h"

#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <tr1/memory>
#include <tr1/unordered_map>
#include <vector>

#include <boost/assert.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

//
// https://en.wikipedia.org/wiki/De_Bruijn_graph
//

template< size_t K >
class DeBruijnGraph {
public:
    typedef std::map< char, size_t > EdgeList;
    struct Node {
        Node() {
        }
        Node(const char edge, size_t weight) {
            children[edge] = weight;
        }
        EdgeList children;
        EdgeList parents;

        size_t indegree() const {
            return parents.size();
        }
        size_t outdegree() const {
            return children.size();
        }

        operator bool() const {
            return (!children.empty()) || (!parents.empty());
        }
    };
    typedef KmerTable< K, Node > NodeList;

    DeBruijnGraph() : _average(0) {
    }

    bool read(std::istream& stream) {
        ARFF arff;

        ARFFReader reader(stream, 2); 
        if (!reader.read(arff)) {
            return false;
        }


        _average = boost::lexical_cast< double >(arff.attrs["mean"]);

        LOG4CXX_DEBUG(logger, boost::format("attrs=[%d],data=[%d] average=[%f]") % arff.attrs.size() % arff.data.size() % _average);

        LOG4CXX_DEBUG(logger, "build debruijn graph begin");
        size_t L = Kmer< K >::length();
        BOOST_FOREACH(const ARFF::ItemPtr& item, arff.data) {
            Kmer< K > prefix((*item)[0], 0, L), suffix((*item)[0], 1, L + 1);

            addEdge(prefix, suffix, boost::lexical_cast< size_t >((*item)[1]));
        }
        LOG4CXX_DEBUG(logger, "build debruijn graph end");

        return true;
    }

    void addEdge(const Kmer< K >& from, const Kmer< K >& to, size_t weight = 1) {
        size_t L = Kmer< K >::length();

        //LOG4CXX_TRACE(logger, boost::format("addEdge: from=[%s], to=[%s], weight=%d") % from % to % weight);

        // children
        {
            typename NodeList::iterator it = _nodelist.find(from);
            char nucleotide = to.nucleotide(L - 1);
            if (it != _nodelist.end()) {
                it->second.children[nucleotide] += weight;
            } else {
                _nodelist.insert(std::make_pair(from, Node(nucleotide, weight)));
            }
        }
        // parents
        {
            typename NodeList::iterator it = _nodelist.find(to);
            char nucleotide = from.nucleotide(0);
            if (it != _nodelist.end()) {
                it->second.parents[nucleotide] += weight;
            } else {
                // Add a virtual node
                Node n;
                n.parents[nucleotide] = weight;
                _nodelist[to] = n;
            }
        }
    }
    void removeEdge(const Kmer< K >& from, const Kmer< K >& to) {
        size_t L = Kmer< K >::length();

        // children
        {
            typename NodeList::iterator i = _nodelist.find(from);
            if (i != _nodelist.end()) {
                char nucleotide = to.nucleotide(L - 1);
                EdgeList::iterator k = i->second.children.find(nucleotide);
                if (k != i->second.children.end()) {
                    i->second.children.erase(k);
                }
            }
        }
        // parents
        {
            typename NodeList::iterator j = _nodelist.find(to);
            if (j != _nodelist.end()) {
                char nucleotide = from.nucleotide(0);
                EdgeList::iterator k = j->second.parents.find(nucleotide);
                if (k != j->second.parents.end()) {
                    j->second.parents.erase(k);
                }
            }
        }
    }
    // Remove suspicious nodes
    void removeNoise();
    double average() const {
        return _average;
    }
private:
    template< size_t X >
    friend std::ostream& operator << (std::ostream& os, const DeBruijnGraph< X >& graph);
    template< size_t X >
    friend class NodeRemover;

    NodeList _nodelist;
    double _average;

    static log4cxx::LoggerPtr logger;
};

template< size_t K >
log4cxx::LoggerPtr DeBruijnGraph< K >::logger(log4cxx::Logger::getLogger("arcs.DeBruijnGraph"));

template< size_t K >
struct KmerCoveragePlus {
    KmerCoveragePlus(const typename DeBruijnGraph< K >::NodeList& nodelist) : _nodelist(nodelist) {
    } 

    double operator()(double l, const Kmer< K >& r) const {
        typename DeBruijnGraph< K >::NodeList::const_iterator i = _nodelist.find(r);
        if (i != _nodelist.end()) {
            for (typename DeBruijnGraph< K >::EdgeList::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
                l += j->second;
            }
        }
        return l;
    }
private:
    const typename DeBruijnGraph< K >::NodeList& _nodelist;
};

//
// Make each node to be an isolated node
//
template< size_t K >
struct NodeRemover {
    NodeRemover(DeBruijnGraph< K >* graph) : _graph(graph) {
    } 

    void add(const Kmer< K >& node) {
        //log4cxx::LoggerPtr logger = DeBruijnGraph< K >::logger;

        //LOG4CXX_TRACE(logger, boost::format("NodeRemover: add(%s)") % node);
        _noiselist.insert(node);
    }
    void add(const Kmer< K >& from, const Kmer< K >& to) {
        //log4cxx::LoggerPtr logger = DeBruijnGraph< K >::logger;

        //LOG4CXX_TRACE(logger, boost::format("NodeRemover: add(%s,%s)") % from % to);
        _edgelist.push_back(std::make_pair(from, to));
    }
    void remove() {
        BOOST_FOREACH(const Kmer< K >& node, _noiselist) {
            unlink(node);
        }
        BOOST_FOREACH(const NoiseEdge& edge, _edgelist) {
            unlink(edge);
        }
    }
    size_t size() const {
        return _noiselist.size();
    }

private:
    typedef std::set< Kmer< K > > NoiseNodeList;
    typedef std::pair< Kmer< K >, Kmer< K > > NoiseEdge;
    typedef std::vector< NoiseEdge > NoiseEdgeList;

    void unlink(const Kmer< K >& kmer) {
        typename DeBruijnGraph< K >::NodeList::iterator k = _graph->_nodelist.find(kmer);
        if (k != _graph->_nodelist.end()) {
            size_t L = Kmer< K >::length();

            // Remove child's links
            {
                int nucleotide = kmer.nucleotide(0);
                for (typename DeBruijnGraph< K >::EdgeList::const_iterator i = k->second.children.begin(); i != k->second.children.end(); ++i) {
                    typename DeBruijnGraph< K >::NodeList::iterator j = _graph->_nodelist.find(kmer + i->first);
                    if (j != _graph->_nodelist.end()) {
                        j->second.parents.erase(nucleotide);
                    }
                }
            }
            k->second.children.clear();

            // Remove parent's links
            {
                int nucleotide = kmer.nucleotide(L - 1);
                for (typename DeBruijnGraph< K >::EdgeList::const_iterator i = k->second.parents.begin(); i != k->second.parents.end(); ++i) {
                    typename DeBruijnGraph< K >::NodeList::iterator j = _graph->_nodelist.find(kmer - i->first);
                    if (j != _graph->_nodelist.end()) {
                        j->second.children.erase(nucleotide);
                    }
                }
            }
            k->second.parents.clear();

            // remove
            _graph->_nodelist.erase(k);
        }
    }

    void unlink(const NoiseEdge& edge) {
        unlink(edge.first, edge.second);
    }

    void unlink(const Kmer< K >& from, const Kmer< K >& to) {
        size_t L = Kmer< K >::length();

        // Children
        {
            typename DeBruijnGraph< K >::NodeList::iterator i = _graph->_nodelist.find(from);
            if (i != _graph->_nodelist.end()) {
                int nucleotide = to.nucleotide(L - 1);
                typename DeBruijnGraph< K >::EdgeList::iterator k = i->second.children.find(nucleotide);
                if (k != i->second.children.end()) {
                    i->second.children.erase(k);
                }
            }
        }
        // Parent
        {
            typename DeBruijnGraph< K >::NodeList::iterator j = _graph->_nodelist.find(to);
            if (j != _graph->_nodelist.end()) {
                int nucleotide = from.nucleotide(0);
                typename DeBruijnGraph< K >::EdgeList::iterator k = j->second.parents.find(nucleotide);
                if (k != j->second.parents.end()) {
                    j->second.parents.erase(k);
                }
            }
        }
        // Remove
        {
            typename DeBruijnGraph< K >::NodeList::iterator i = _graph->_nodelist.find(from);
            if (i != _graph->_nodelist.end() && !i->second) {
                _graph->_nodelist.erase(i);
            }
            typename DeBruijnGraph< K >::NodeList::iterator j = _graph->_nodelist.find(to);
            if (j != _graph->_nodelist.end() && !j->second) {
                _graph->_nodelist.erase(j);
            }
        }
    }

    NoiseNodeList _noiselist;
    NoiseEdgeList _edgelist;

    DeBruijnGraph< K >* _graph;
};


template< size_t K >
void DeBruijnGraph< K >::removeNoise() {
    LOG4CXX_DEBUG(logger, boost::format("remove noise begin. nodelist=[%d]") % _nodelist.size());

    NodeRemover< K > remover(this);

    for (typename NodeList::const_iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
        // Condensed node
        if ((i->second.indegree() == 1 && i->second.outdegree() == 1)) {
            continue;
        }

        for (typename EdgeList::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
            std::list< Kmer< K > > group; // h-path

            // Look forward
            for (typename NodeList::const_iterator k = _nodelist.find(i->first + j->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != i->first; k = _nodelist.find(k->first + k->second.children.begin()->first)) {
                group.push_back(k->first);
            }

            // front && back
            typename NodeList::const_iterator front(i), back(group.empty() ? _nodelist.find(i->first + j->first) : _nodelist.find(group.back() + _nodelist.find(group.back())->second.children.begin()->first));
            BOOST_ASSERT(front != _nodelist.end() && back != _nodelist.end());

            // length && coverage
            size_t length = 1 + group.size() + 1; // first + group + last
            double coverage = ((double)(j->second + std::accumulate(group.begin(), group.end(), 0, KmerCoveragePlus< K >(_nodelist)))) / (length - 1); 

            //LOG4CXX_TRACE(logger, boost::format("removeNoise: front=[%s],length=[%d],coverage=[%f],group=[%d],front=[indegree(%d),outdegree(%d)],back=[indegree(%d),outdegree(%d)]") % front->first % length % coverage % group.size() % front->second.indegree() % front->second.outdegree() % back->second.indegree() % back->second.outdegree());

            // rules
            if ((front->second.indegree() == 0 && front->second.outdegree() == 1) || (back->second.indegree() == 1 && back->second.outdegree() == 0)) { // tips
                if (K + length < 2 * (K + 1) || coverage < _average / 5.0) {
                    if (front->second.indegree() == 0 && front->second.outdegree() == 1) {
                        remover.add(front->first);
                    }
                    if (back->second.indegree() == 1 && back->second.outdegree() == 0) {
                        remover.add(back->first);
                    }
                    std::for_each(group.begin(), group.end(), boost::bind(&NodeRemover< K >::add, &remover, _1));
                    if (group.empty()) {
                        remover.add(front->first, back->first);
                    }
                }
            } else if ((K + length < 2 * (K + 1) && coverage < _average / 5.0) || (coverage <= 3.0) || (coverage < _average / 10.0)) { // buble or link ?????????
                std::for_each(group.begin(), group.end(), boost::bind(&NodeRemover< K >::add, &remover, _1));
                if (group.empty()) {
                    remover.add(front->first, back->first);
                }
            }
        }
    }

    remover.remove();

    LOG4CXX_DEBUG(logger, boost::format("remove %d noise nodes, %d nodes left") % remover.size() % _nodelist.size());
}

template< size_t K >
struct KmerIndexer {
    typedef KmerTable< K, size_t > IndexTable;

    KmerIndexer(const typename DeBruijnGraph< K >::NodeList& nodelist) {
        _index_tbl.rehash(nodelist.size());

        size_t i = 0;
        for (typename DeBruijnGraph< K >::NodeList::const_iterator it = nodelist.begin(); it != nodelist.end(); ++it) {
            _index_tbl[it->first] = ++i;
        }
    }
    size_t operator[](const Kmer< K >& kmer) const {
        typename IndexTable::const_iterator it = _index_tbl.find(kmer);
        if (it != _index_tbl.end()) {
            return it->second;
        }
        return 0;
    }
private:
    IndexTable _index_tbl;
};

template< size_t K >
std::ostream& operator << (std::ostream& os, const DeBruijnGraph< K >& graph) {
    KmerIndexer< K > indexer(graph._nodelist);

    size_t L = Kmer< K >::length();
    // Merge nodes with indegree == outdegree == 1
    for (typename DeBruijnGraph< K >::NodeList::const_iterator i = graph._nodelist.begin(); i != graph._nodelist.end(); ++i) {
        // edge
        if (i->second.indegree() == 1 && i->second.outdegree() == 1) {
            continue;
        }
        for (typename DeBruijnGraph< K >::EdgeList::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
            std::string edge = i->first.sequence() + j->first;
            size_t length = 1; // first
            double coverage = j->second;

            typename DeBruijnGraph< K >::NodeList::const_iterator k = graph._nodelist.find(i->first + j->first);
            while (k != graph._nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1) {
                length += 1;
                coverage += k->second;
                k = graph._nodelist.find(k->first + k->second.children.begin()->first);
                if (k != graph._nodelist.end()) {
                    edge += k->first.nucleotide(L - 1);
                }
            }
            os << boost::format("%d\t%d\t%s\n") % indexer[i->first] % indexer[i->first + j->first] % edge;
            os << boost::format("%f\n") % (coverage / length);
        }
    }

    return os;
}

#endif // debruijn_graph_h_
