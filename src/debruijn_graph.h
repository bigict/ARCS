#ifndef debruijn_graph_h_
#define debruijn_graph_h_

#include "arff.h"
#include "kmer.h"

#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
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
    typedef std::unordered_map< char, uint32_t > EdgeList;
    struct Node {
        Node() : children(Nucleotide::NUM), parents(Nucleotide::NUM) {
        }
        Node(const char edge, size_t weight) : children(Nucleotide::NUM), parents(Nucleotide::NUM) {
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

    bool read(std::istream& stream, size_t iter) {
        ARFF arff;

        ARFFReader reader(stream, 2); 
        if (!reader.read(arff)) {
            return false;
        }


        _average = boost::lexical_cast< double >(arff.attrs["mean"]);

        LOG4CXX_DEBUG(logger, boost::format("attrs=[%d],data=[%d] average=[%f]") % arff.attrs.size() % arff.data.size() % _average);

        _nodelist.rehash(arff.data.size() * 2);

        LOG4CXX_DEBUG(logger, "build debruijn graph begin");
        size_t L = Kmer< K >::length();
        size_t usedKmerNum = 0;
        BOOST_FOREACH(const ARFF::ItemPtr& item, arff.data) {
            Kmer< K > prefix((*item)[0], 0, L), suffix((*item)[0], 1, L + 1);

            size_t copyNum = boost::lexical_cast< size_t >((*item)[1]);
            if(iter == 1) {
                if(copyNum != 1) {
                    usedKmerNum += addEdge(prefix, suffix, copyNum); 
                }
            } else if(copyNum == 1) {
                usedKmerNum += addEdge(prefix, suffix, copyNum); 
            }
        }
        LOG4CXX_DEBUG(logger, boost::format("used kmer number = %d") % usedKmerNum);
        LOG4CXX_DEBUG(logger, "build debruijn graph end");

        return true;
    }

    int addEdge(const Kmer< K >& from, const Kmer< K >& to, size_t weight = 1) {
        size_t L = Kmer< K >::length();

        //LOG4CXX_TRACE(logger, boost::format("addEdge: from=[%s], to=[%s], weight=%d") % from % to % weight);
        if(weight == 1) {
            typename NodeList::iterator it = _nodelist.find(from);
            typename NodeList::iterator jt = _nodelist.find(to);
            /*
            if(it != _nodelist.end() && jt != _nodelist.end()
                && (it->second.children.size() != 0 || jt->second.parents.size() != 0)) {
                return 0;
            }
            */
            if(it == _nodelist.end() || jt == _nodelist.end() || 
                    it->second.children.size() != 0 || jt->second.parents.size() != 0) { // only when outdegree of it == 0 and indegree of jt == 0, add edge
                return 0;
            }
        } 
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
        return 1;
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
    LOG4CXX_DEBUG(logger, boost::format("remove noise begin. nodelist=[%d] _K=[%d]") % _nodelist.size() % Kmer< K >::length());

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
                if (Kmer< K >::length() + length < 2 * (Kmer< K >::length() + 1) || coverage < _average / 5.0) {
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
            } else if ((Kmer< K >::length() + length < 2 * (Kmer< K >::length() + 1) && coverage < _average / 5.0) || (coverage <= 3.0) || (coverage < _average / 10.0)) { // buble or link ?????????
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
    std::vector< size_t > contigLength;
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
                coverage += k->second.children.begin()->second;
                k = graph._nodelist.find(k->first + k->second.children.begin()->first);
                if (k != graph._nodelist.end()) {
                    edge += k->first.nucleotide(L - 1);
                }
            }
            os << boost::format("%d\t%d\t%s\n") % indexer[i->first] % indexer[i->first + j->first] % edge;
            os << boost::format("%f\n") % (coverage / length);
            contigLength.push_back( edge.length() );
        }
    }
    std::sort(contigLength.begin(), contigLength.end(), std::greater<int>());
    size_t sumLength = std::accumulate(contigLength.begin(), contigLength.end(), 0);
    size_t N50 = sumLength * 0.5;
    size_t N90 = sumLength * 0.9;
    size_t sum = 0, N50_val, N90_val;
    bool findN50 = false;
    for(std::vector< size_t >::iterator it = contigLength.begin(); it != contigLength.end(); ++it) {
        sum += (*it);
        if(!findN50 && sum >= N50) {
            N50_val = (*it);
            findN50 = true;
        }
        if(sum >= N90) {
            N90_val = (*it);
            break;
        }
    }
	LOG4CXX_INFO(DeBruijnGraph< K >::logger, boost::format("contig N50 : %d") % N50_val);
    LOG4CXX_INFO(DeBruijnGraph< K >::logger, boost::format("contig N90 : %d") % N90_val);
    LOG4CXX_INFO(DeBruijnGraph< K >::logger, boost::format("contig number : %d") % contigLength.size());
    
    return os;
}

#endif // debruijn_graph_h_
