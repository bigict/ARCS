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
    typedef std::map< size_t, size_t > EdgeList;
    struct Node {
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
    typedef std::vector< Node > NodeList;

    DeBruijnGraph() : _average(0) {
    }

    bool read(std::istream& stream) {
        ARFF arff;

        ARFFReader reader(stream); 
        if (!reader.read(arff)) {
            return false;
        }

        LOG4CXX_DEBUG(logger, boost::format("attrs=[%d],data=[%d]") % arff.attrs.size() % arff.data.size());

        _indexer.init(arff);

        _average = boost::lexical_cast< double >(arff.attrs["mean"]);
        _nodelist.resize(_indexer.count());

        LOG4CXX_DEBUG(logger, boost::format("average=[%f],nodes=[%d]") % _average % _nodelist.size());

        LOG4CXX_DEBUG(logger, "build debruijn graph begin");
        BOOST_FOREACH(const ARFF::Item& item, arff.data) {
            Kmer< K > prefix(item[0], 0, K), suffix(item[0], 1, K + 1);

            addEdge(_indexer.id(prefix), _indexer.id(suffix), boost::lexical_cast< size_t >(item[1]));
        }
        LOG4CXX_DEBUG(logger, "build debruijn graph end");

        return true;
    }

    void addEdge(size_t i, size_t j, size_t weight = 1) {
        Node& from(_nodelist[i]);
        Node& to(_nodelist[j]);
        // children
        {
            EdgeList::iterator k = from.children.find(j);
            if (k != from.children.end()) {
                k->second += weight;
            } else {
                from.children[j] = weight;
            }
        }
        // parents
        {
            EdgeList::iterator k = to.parents.find(i);
            if (k != to.parents.end()) {
                k->second += weight;
            } else {
                to.parents[i] = weight;
            }
        }
    }
    void removeEdge(size_t i, size_t j) {
        Node& from(_nodelist[i]);
        Node& to(_nodelist[j]);
        // children
        {
            EdgeList::iterator k = from.children.find(j);
            if (k != from.children.end()) {
                from.children.erase(k);
            }
        }
        // parents
        {
            EdgeList::iterator k = to.parents.find(i);
            if (k != to.parents.end()) {
                to.parents.erase(k);
            }
        }
    }
    // Remove suspicious nodes
    void removeNoise();
private:
    template< size_t X >
    friend std::ostream& operator << (std::ostream& os, const DeBruijnGraph< X >& graph);

    class Indexer {
    public:
        typedef std::vector< Kmer< K > > KmerList;

        Indexer() {
        }
        Indexer(const ARFF& arff) {
            init(arff);
        }

        void init(const ARFF& arff) {
            _data.clear();

            BOOST_FOREACH(const ARFF::Item& item, arff.data) {
                BOOST_ASSERT(item.size() == 2);
                BOOST_ASSERT(item[0].length() == K + 1);

                Kmer< K > prefix(item[0], 0, K), suffix(item[0], 1, K + 1);
                _data.push_back(prefix);
                _data.push_back(suffix);

            }

            // uniq
            std::sort(_data.begin(), _data.end());
            typename KmerList::iterator last = std::unique(_data.begin(), _data.end());
            _data.resize(std::distance(_data.begin(), last));
        }

        size_t id(const Kmer< K >& kmer) const {
            typename KmerList::const_iterator i = std::lower_bound(_data.begin(), _data.end(), kmer);
            if (*i == kmer) {
                return std::distance(_data.begin(), i);
            }
            return -1;
        }
        const Kmer< K >& kmer(size_t id) const {
            return _data[id];
        }

        size_t count() const {
            return _data.size();
        }

    private:
        KmerList _data;
    };

    NodeList _nodelist;
    Indexer _indexer;

    double _average;

    static log4cxx::LoggerPtr logger;
};

template< size_t K >
log4cxx::LoggerPtr DeBruijnGraph< K >::logger(log4cxx::Logger::getLogger("arcs.DeBruijnGraph"));

template< size_t K >
struct KmerCoveragePlus {
    KmerCoveragePlus(const typename DeBruijnGraph< K >::NodeList& nodelist) : _nodelist(nodelist) {
    } 

    double operator()(double l, size_t r) const {
        const typename DeBruijnGraph< K >::Node& node(_nodelist[r]);
        for (typename DeBruijnGraph< K >::EdgeList::const_iterator j = node.children.begin(); j != node.children.end(); ++j) {
            l += j->second;
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
    NodeRemover(typename DeBruijnGraph< K >::NodeList& nodelist) : _nodelist(nodelist) {
    } 

    void add(size_t node) {
        _noiselist.insert(node);
    }
    void add(size_t i, size_t j) {
        _edgelist.push_back(std::make_pair(i, j));
    }
    void remove() {
        BOOST_FOREACH(size_t node, _noiselist) {
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
    typedef std::set< size_t > NoiseNodeList;
    typedef std::pair< size_t, size_t > NoiseEdge;
    typedef std::vector< NoiseEdge > NoiseEdgeList;

    void unlink(size_t n) {
        typename DeBruijnGraph< K >::Node& node = _nodelist[n];

        // Remove child's links
        for (typename DeBruijnGraph< K >::EdgeList::const_iterator i = node.children.begin(); i != node.children.end(); ++i) {
            typename DeBruijnGraph< K >::Node& j = _nodelist[i->first];
            typename DeBruijnGraph< K >::EdgeList::iterator k = j.parents.find(n);
            if (k != j.parents.end()) {
                j.parents.erase(n);
            }
        }
        node.children.clear();

        // Remove parent's links
        for (typename DeBruijnGraph< K >::EdgeList::const_iterator i = node.parents.begin(); i != node.parents.end(); ++i) {
            typename DeBruijnGraph< K >::Node& j = _nodelist[i->first];
            typename DeBruijnGraph< K >::EdgeList::iterator k = j.children.find(n);
            if (k != j.parents.end()) {
                j.children.erase(k);
            }
        }
        node.parents.clear();
    }

    void unlink(const NoiseEdge& edge) {
        unlink(edge.first, edge.second);
    }

    void unlink(size_t i, size_t j) {
        // Children
        {
            typename DeBruijnGraph< K >::Node& from = _nodelist[i];
            typename DeBruijnGraph< K >::EdgeList::iterator k = from.children.find(j);
            if (k != from.children.end()) {
                from.children.erase(k);
            }
        }
        // Parent
        {
            typename DeBruijnGraph< K >::Node& to = _nodelist[j];
            typename DeBruijnGraph< K >::EdgeList::iterator k = to.parents.find(i);
            if (k != to.parents.end()) {
                to.parents.erase(k);
            }
        }
    }

    NoiseNodeList _noiselist;
    NoiseEdgeList _edgelist;

    typename DeBruijnGraph< K >::NodeList& _nodelist;
};


template< size_t K >
void DeBruijnGraph< K >::removeNoise() {
    LOG4CXX_DEBUG(logger, boost::format("remove noise begin. nodelist=[%d]") % _nodelist.size());

    NodeRemover< K > remover(_nodelist);

    for (size_t i = 0; i < _nodelist.size(); ++i) {
        Node& node = _nodelist[i];

        // Condensed node
        if ((node.indegree() == 1 && node.outdegree() == 1)) {
            continue;
        }

        for (EdgeList::const_iterator j = node.children.begin(); j != node.children.end(); ++j) {
            std::list< size_t > group; // h-path

            // Look forward
            for (size_t k = j->first; _nodelist[k].indegree() == 1 && _nodelist[k].outdegree() == 1 && k != i; k = _nodelist[k].children.begin()->first) {
                group.push_back(k);
            }

            // front && back
            size_t front(i), back(group.empty() ? j->first : _nodelist[group.back()].children.begin()->first);

            // length && coverage
            size_t length = 1 + group.size() + 1; // first + group + last
            double coverage = ((double)(j->second + std::accumulate(group.begin(), group.end(), 0, KmerCoveragePlus< K >(_nodelist)))) / (length - 1); 

            LOG4CXX_TRACE(logger, boost::format("removeNoise: front=[%s],length=[%d],coverage=[%f],group=[%d],front=[indegree(%d),outdegree(%d)],back=[indegree(%d),outdegree(%d)]") % _indexer.kmer(front) % length % coverage % group.size() % _nodelist[front].indegree() % _nodelist[front].outdegree() % _nodelist[back].indegree() % _nodelist[back].outdegree());

            // rules
            if ((_nodelist[front].indegree() == 0 && _nodelist[front].outdegree() == 1) || (_nodelist[back].indegree() == 1 && _nodelist[back].outdegree() == 0)) { // tips
                if (K + length < 2 * (K + 1) || coverage < _average / 5.0) {
                    if (_nodelist[front].indegree() == 0 && _nodelist[front].outdegree() == 1) {
                        remover.add(front);
                    }
                    if (_nodelist[back].indegree() == 1 && _nodelist[back].outdegree() == 0) {
                        remover.add(back);
                    }
                    std::for_each(group.begin(), group.end(), boost::bind(&NodeRemover< K >::add, &remover, _1));
                    if (group.empty()) {
                        remover.add(front, back);
                    }
                }
            } else if ((K + length < 2 * (K + 1) && coverage < _average / 5.0) || (coverage <= 3.0) || (coverage < _average / 10.0)) { // buble or link ?????????
                std::for_each(group.begin(), group.end(), boost::bind(&NodeRemover< K >::add, &remover, _1));
                if (group.empty()) {
                    remover.add(front, back);
                }
            }
        }
    }

    remover.remove();

    LOG4CXX_DEBUG(logger, boost::format("remove noise end. nodelist=[%d-%d=%d]") % _nodelist.size() % remover.size() % (_nodelist.size() - remover.size()));
}

template< size_t K >
std::ostream& operator << (std::ostream& os, const DeBruijnGraph< K >& graph) {
    // Merge nodes with indegree == outdegree == 1
    for (size_t i = 0; i < graph._nodelist.size(); ++i) {
        const typename DeBruijnGraph< K >::Node& node = graph._nodelist[i];

        // edge
        if (node.indegree() == 1 && node.outdegree() == 1) {
            continue;
        }
        for (typename DeBruijnGraph< K >::EdgeList::const_iterator j = node.children.begin(); j != node.children.end(); ++j) {
            std::string edge = graph._indexer.kmer(i).sequence() + graph._indexer.kmer(j->first).nucleotide(K - 1);
            size_t length = 1; // first
            double coverage = j->second;

            size_t k = j->first;
            while (graph._nodelist[k].indegree() == 1 && graph._nodelist[k].outdegree() == 1) {
                const typename DeBruijnGraph< K >::EdgeList::const_iterator x = graph._nodelist[k].children.begin();

                length += 1;
                coverage += x->second;
                edge += graph._indexer.kmer(x->first).nucleotide(K - 1);

                k = x->first;
            }
            os << boost::format("%d\t%d\t%s\t%f\t%d") % i % k % edge % coverage % length << std::endl;
            os << (coverage / length) << std::endl;
        }
    }

    return os;
}

#endif // debruijn_graph_h_
