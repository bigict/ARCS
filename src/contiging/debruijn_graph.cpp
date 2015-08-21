#include "debruijn_graph.h"
#include "kmer_tbl.h"

#include <algorithm>
#include <list>
#include <numeric>
#include <set>
#include <fstream>
#include <tr1/unordered_set>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.debruijn"));

DeBruijnGraph::DeBruijnGraph(const KmerTable& tbl) : _K(tbl.K()), _average(0) {
    LOG4CXX_DEBUG(logger, boost::format("build debruijn graph from kmer hash table with K=[%d]") % _K);

    double delta;
    tbl.statistics(&_average, &delta);

    LOG4CXX_DEBUG(logger, boost::format("average coverage=[%f],delta=[%f]") % _average % delta);

    tbl.buildDeBruijn(this);
}

DeBruijnGraph::~DeBruijnGraph() {
}

void DeBruijnGraph::addKmer(const Kmer& kmer, size_t weight) {
    size_t k = kmer.length();
    Kmer node = kmer.subKmer(0, k - 1);
    Kmer edge = kmer.subKmer(1);

    LOG4CXX_TRACE(logger, boost::format("addKmer: %s with node=[%s], edge=[%s]") % kmer % node % edge);

    // children
    {
        NodeList::iterator it = _nodelist.find(node);
        if (it != _nodelist.end()) {
            it->second.children[edge] += weight;
        } else {
            _nodelist[node] = Node(edge, weight);
        }
    }
    // parents
    {
        NodeList::iterator it = _nodelist.find(edge);
        if (it != _nodelist.end()) {
            it->second.parents[node] += weight;
        } else {
            // Add a virtual node
            Node n;
            n.parents[node] = weight;
            _nodelist[edge] = n;
        }
    }
}

void DeBruijnGraph::removeKmer(const Kmer& kmer) {
    size_t k = kmer.length();
    Kmer node = kmer.subKmer(0, k - 1);
    Kmer edge = kmer.subKmer(1);

    LOG4CXX_TRACE(logger, boost::format("removeKmer: %s with key=[%s], edge=[%s]") % kmer % node % edge);
    
    // TODO: bugfix
    NodeList::iterator it = _nodelist.find(node);
    if (it != _nodelist.end()) {
        it->second.children.erase(edge);
        if (!it->second) {
            _nodelist.erase(it);
        }
    }
}

void DeBruijnGraph::removeEdge(const Kmer& src, const Kmer& dest) {
    LOG4CXX_TRACE(logger, boost::format("removeEdge: %s src=[%s], dest=[%s]") % src % dest);

    NodeList::iterator i = _nodelist.find(src);
    if (i != _nodelist.end()) {
        EdgeList::iterator j = i->second.children.find(dest);
        if (j != i->second.children.end()) {
            i->second.children.erase(j);
        }
    }
}

struct KmerLengthPlus {
    KmerLengthPlus(size_t K) : _K(K) {
    }
    size_t operator()(size_t l, const Kmer& r) const {
        return l + length(_K, r);
    }

    static size_t length(size_t K, const Kmer& kmer) {
        return kmer.length() - (K-1) + 1;
    }
private:
    size_t _K;
};

struct KmerCoveragePlus {
    KmerCoveragePlus(const DeBruijnGraph::NodeList* nodelist, size_t K) : _nodelist(nodelist), _K(K) {
    } 

    double operator()(double l, const Kmer& r) const {
        DeBruijnGraph::NodeList::const_iterator i = _nodelist->find(r);
        if (i != _nodelist->end()) {
            for (DeBruijnGraph::EdgeList::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
                l += coverage(_K, r, j->second);
            }
        }
        return l;
    }

    static double coverage(size_t K, const Kmer& kmer, double n) {
        return n * KmerLengthPlus::length(K, kmer);
    }
private:
    const DeBruijnGraph::NodeList* _nodelist;
    size_t _K;
};

void DeBruijnGraph::compact() {
    //typedef std::tr1::unordered_set< Kmer, KmerHasher > CountingList;
    typedef std::set< Kmer > CountingList;

    LOG4CXX_DEBUG(logger, boost::format("compact with %d nodes, begin") % _nodelist.size());

    // Searching nodes with indegree==1 and out outdegree==1
    CountingList merge_nodelist;
    {
        LOG4CXX_DEBUG(logger, boost::format("counting indegree and outdegree begin"));

        // Counting in-degrees & out-degrees
        for (NodeList::const_iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
            LOG4CXX_TRACE(logger, boost::format("kmer=[%s],indegree=[%d],outdegree=[%d]") % i->first % i->second.indegree() % i->second.outdegree());
            if (i->second.indegree() == 1 && i->second.outdegree() == 1) {
                merge_nodelist.insert(i->first);
            }
        }

        LOG4CXX_DEBUG(logger, boost::format("counting indegree and outdegree end"));
    }

    LOG4CXX_DEBUG(logger, boost::format("%d nodes found") % merge_nodelist.size());

    // Merge
    while (!merge_nodelist.empty()) {
        CountingList::iterator i = merge_nodelist.begin();
        
        NodeList::const_iterator j = _nodelist.find(*i);
        BOOST_ASSERT(j != _nodelist.end());

        std::list< Kmer > group = boost::assign::list_of(j->first);

        // Look backword.
        for (NodeList::const_iterator k = _nodelist.find(j->second.parents.begin()->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != j->first && k->second.parents.find(j->first) == k->second.parents.end(); k = _nodelist.find(k->second.parents.begin()->first)) {
            group.push_front(k->first);
        }

        // Look forward
        for (NodeList::const_iterator k = _nodelist.find(j->second.children.begin()->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != j->first && k->second.children.find(j->first) == k->second.children.end(); k = _nodelist.find(k->second.children.begin()->first)) {
            group.push_back(k->first);
        }

        // Merge nodes
        if (group.size() > 1) {
            Kmer key;

            BOOST_FOREACH(const Kmer& kmer, group) {
                if (key.length() == 0) {
                    key = kmer;
                } else {
                    key += kmer.subKmer(_K - 2);//key->kmer????
                }
                LOG4CXX_TRACE(logger, boost::format("compact: key=[%s]/[%s]") % key % kmer);
            }

            Node val;

            // Average coverage
            size_t length = std::accumulate(group.begin(), group.end(), 0, KmerLengthPlus(_K));
            double coverage = std::accumulate(group.begin(), group.end(), 0, KmerCoveragePlus(&_nodelist, _K));

            // Set parrents
            Node& front(_nodelist[group.front()]);
            val.parents = front.parents;
            for (EdgeList::const_iterator it = front.parents.begin(); it != front.parents.end(); ++it) {
                NodeList::iterator k = _nodelist.find(it->first);
                if (k != _nodelist.end()) {
                    k->second.children[key] = k->second.children[group.front()];
                    k->second.children.erase(group.front());
                }
            }

            // Set children
            Node& back(_nodelist[group.back()]);
            val.children = back.children;
            for (EdgeList::const_iterator it = back.children.begin(); it != back.children.end(); ++it) {

                val.children[it->first] = coverage / length; // Use average coverage here.

                NodeList::iterator k = _nodelist.find(it->first);
                if (k != _nodelist.end()) {
                    k->second.parents[key] = k->second.parents[group.back()];
                    k->second.parents.erase(group.back());
                }
            }

            BOOST_FOREACH(const Kmer& kmer, group) {
                _nodelist.erase(kmer);
            }
            
            // Add a new condensed node
            _nodelist[key] = val;
        }
        BOOST_FOREACH(const Kmer& kmer, group) {
            merge_nodelist.erase(kmer);
        }

        LOG4CXX_TRACE(logger, boost::format("merge %d nodes with %d nodes left to be merged")  % group.size()% merge_nodelist.size());
    }

    LOG4CXX_DEBUG(logger, boost::format("compact with %d nodes, end") % _nodelist.size());
}

//
// Make each node to be an isolated node
//
struct KmerRemover {
    KmerRemover(DeBruijnGraph::NodeList* nodelist) : _nodelist(nodelist) {
    } 

    void transaction_begin() {
        //LOG4CXX_TRACE(logger, boost::format("KmerRemover: transaction_begin"));
        //_translist.clear();
    }
    void transaction_add(const Kmer& kmer) {
        LOG4CXX_DEBUG(logger, boost::format("KmerRemover: transaction_add, kmer=[%s]") % kmer);
        _noiselist.push_back(kmer);
        //_translist.push_back(kmer);
    }
    void transaction_add_edge(const Kmer& from, const Kmer& to) {
        _edgelist.push_back(std::make_pair(from, to));
        //LOG4CXX_DEBUG(logger, boost::format("KmerRemover: transaction_add_edge, from=[%s],to=[%s]/edgelist=[%d]") % from % to % _edgelist.size());
    }
    void transaction_end() {
        //std::for_each(_translist.begin(), _translist.end(), boost::bind(&KmerRemover::unlink, this, _1));
        //LOG4CXX_TRACE(logger, boost::format("KmerRemover: transaction_end, nodes=[%d]") % _translist.size());
    }
    
    void remove() {
        /*
        BOOST_FOREACH(const Kmer& kmer, _noiselist) {
            unlink(kmer);
        }
        */
        std::for_each(_noiselist.begin(), _noiselist.end(), boost::bind(&KmerRemover::unlink, this, _1));
        for (NoiseEdgeList::const_iterator i = _edgelist.begin(); i != _edgelist.end(); ++i) {
            unlink(i->first, i->second);
        }
    }
    size_t size() const {
        return _noiselist.size();
    }

private:
    void unlink(const Kmer& kmer) {
        DeBruijnGraph::NodeList::iterator k = _nodelist->find(kmer);
        if (k != _nodelist->end()) {
            // Remove child's links
            for (DeBruijnGraph::EdgeList::const_iterator i = k->second.children.begin(); i != k->second.children.end(); ++i) {
                DeBruijnGraph::NodeList::iterator j = _nodelist->find(i->first);
                if (j != _nodelist->end()) {
                    j->second.parents.erase(kmer);
                }
            }
            k->second.children.clear();

            // Remove parent's links
            for (DeBruijnGraph::EdgeList::const_iterator i = k->second.parents.begin(); i != k->second.parents.end(); ++i) {
                DeBruijnGraph::NodeList::iterator j = _nodelist->find(i->first);
                if (j != _nodelist->end()) {
                    j->second.children.erase(kmer);
                }
            }
            k->second.parents.clear();
            
            // Add it to noise node list
            // _noiselist.push_back(kmer);
        }
        _nodelist->erase(kmer);
    }

    void unlink(const Kmer& from, const Kmer& to) {
        LOG4CXX_DEBUG(logger, boost::format("Remove: remove edge(%s,%s)") % from % to);
        // Children
        {
            DeBruijnGraph::NodeList::iterator i = _nodelist->find(from);
            if (i != _nodelist->end()) {
                DeBruijnGraph::EdgeList::iterator k = i->second.children.find(to);
                if (k != i->second.children.end()) {
                    i->second.children.erase(k);
                }
            }
        }
        // Parent
        {
            DeBruijnGraph::NodeList::iterator j = _nodelist->find(to);
            if (j != _nodelist->end()) {
                DeBruijnGraph::EdgeList::iterator k = j->second.parents.find(from);
                if (k != j->second.parents.end()) {
                    j->second.parents.erase(k);
                }
            }
        }
        // Remove
        {
            DeBruijnGraph::NodeList::iterator i = _nodelist->find(from);
            if (i != _nodelist->end() && !i->second) {
                _nodelist->erase(i);
                LOG4CXX_DEBUG(logger, boost::format("Remover: remove kmer=[%s]") % i->first);
            }
            DeBruijnGraph::NodeList::iterator j = _nodelist->find(to);
            if (j != _nodelist->end() && !j->second) {
                _nodelist->erase(j);
                LOG4CXX_DEBUG(logger, boost::format("Remover: remove kmer=[%s]") % j->first);
            }
        }
    }

    typedef std::list< Kmer > NoiseNodeList;
    NoiseNodeList _noiselist;
//    NoiseList _translist;
    DeBruijnGraph::NodeList* _nodelist;
    typedef std::pair< Kmer, Kmer > NoiseEdge;
    typedef std::vector< NoiseEdge > NoiseEdgeList;
    NoiseEdgeList _edgelist;
};


//
// Remove buble and tips. Note that this function is compatible with condensed and uncondensed graph.
//
void DeBruijnGraph::removeNoise() {
    LOG4CXX_DEBUG(logger, boost::format("remove noise with %d nodes, begin") % _nodelist.size());

    KmerRemover remover(&_nodelist);

    for (NodeList::const_iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
        // Condensed node
        if (i->second.indegree() == 1 && i->second.outdegree() == 1) {
            continue;
        }

        EdgeList children = i->second.children; // copy
        for (EdgeList::const_iterator j = children.begin(); j != children.end(); ++j) {
            remover.transaction_begin();

            std::list< Kmer > group;
            // Look forward
            for (NodeList::const_iterator k = _nodelist.find(j->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != i->first; k = _nodelist.find(k->second.children.begin()->first)) {
                group.push_back(k->first);
            }

            // front && back
            NodeList::const_iterator front(i), back(group.empty() ? _nodelist.find(j->first) : _nodelist.find(_nodelist.find(group.back())->second.children.begin()->first));

            // length && coverage
            size_t length = KmerLengthPlus::length(_K, front->first) + std::accumulate(group.begin(), group.end(), 0, KmerLengthPlus(_K)) + KmerLengthPlus::length(_K, back->first);
            double coverage = (KmerCoveragePlus::coverage(_K, i->first, j->second) + std::accumulate(group.begin(), group.end(), 0, KmerCoveragePlus(&_nodelist, _K))) / (length-1); 
            LOG4CXX_DEBUG(logger, boost::format("removeNoise: front=[%s],length=[%d],coverage=[%f],group=[%d],front=[indegree(%d),outdegree(%d)],back=[indegree(%d),outdegree(%d)]") % front->first % length % coverage % group.size() % front->second.indegree() % front->second.outdegree() % back->second.indegree() % back->second.outdegree());

            // rules
            if ((front->second.indegree() == 0 && front->second.outdegree() == 1) || (back->second.indegree() == 1 && back->second.outdegree() == 0)) { // tips
                if (_K + length - 1 < 2 * _K || coverage < _average / 5.0) {
                    if (front->second.indegree() == 0 && front->second.outdegree() == 1) {
                        remover.transaction_add(front->first);
                    }
                    if (back->second.indegree() == 1 && back->second.outdegree() == 0) {
                        remover.transaction_add(back->first);
                    }
                    std::for_each(group.begin(), group.end(), boost::bind(&KmerRemover::transaction_add, &remover, _1));
                    if (group.empty()) {
                        remover.transaction_add_edge(front->first, back->first);
                    }
                }
            } else if ((_K + length - 1 < 2 * _K && coverage < _average / 5.0) || (coverage <= 3.0) || (coverage < _average / 10.0)) { // buble or link ?????????
                std::for_each(group.begin(), group.end(), boost::bind(&KmerRemover::transaction_add, &remover, _1));
                if (group.empty()) {
                    remover.transaction_add_edge(front->first, back->first);
                }
            }
            remover.transaction_end();
        }
    }
    
    remover.remove();

    LOG4CXX_DEBUG(logger, boost::format("remove %d noise nodes, %d nodes left") % remover.size() % _nodelist.size());
}

struct KmerIndexer {
    typedef std::tr1::unordered_map< Kmer, size_t, KmerHasher > IndexTable;

    KmerIndexer(const DeBruijnGraph::NodeList& nodelist) {
        size_t i = 0;
        for (DeBruijnGraph::NodeList::const_iterator it = nodelist.begin(); it != nodelist.end(); ++it) {
            _index_tbl[it->first] = ++i;
        }
    }
    size_t operator[](const Kmer& kmer) const {
        IndexTable::const_iterator it = _index_tbl.find(kmer);
        if (it != _index_tbl.end()) {
            return it->second;
        }
        return 0;
    }
private:
    IndexTable _index_tbl;
};

std::ostream& operator << (std::ostream& os, const DeBruijnGraph& graph) {
    KmerIndexer indexer(graph._nodelist);
    size_t contig_index = 0, component_index = 0;
    
    // Merge nodes with indegree == outdegree == 1
    for (DeBruijnGraph::NodeList::const_iterator i = graph._nodelist.begin(); i != graph._nodelist.end(); ++i) {
        // edge
        if (i->second.indegree() == 1 && i->second.outdegree() == 1) {
            continue;
        }
        for (DeBruijnGraph::EdgeList::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
            Kmer edge = i->first + j->first.subKmer(graph._K - 2);
            size_t length = KmerLengthPlus::length(graph._K, i->first);
            double coverage = KmerCoveragePlus::coverage(graph._K, i->first, j->second);

            std::vector< double > X = boost::assign::list_of(j->second);

            DeBruijnGraph::NodeList::const_iterator k = graph._nodelist.find(j->first);
            while (k != graph._nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1) {
                length += KmerLengthPlus::length(graph._K, k->first);
                coverage += KmerCoveragePlus::coverage(graph._K, k->first, k->second.children.begin()->second);
                X.push_back(KmerCoveragePlus::coverage(graph._K, k->first, k->second.children.begin()->second));
                k = graph._nodelist.find(k->second.children.begin()->first);
                if (k != graph._nodelist.end()) {
                    edge += k->first.subKmer(graph._K - 2);
                }
            }

            os << boost::format("%d\t%d\t%s\t%f\t%f\t%d\t%d") % indexer[i->first] % indexer[j->first] % edge % (coverage / length) % coverage % length % i->second.indegree() << std::endl;
            BOOST_FOREACH(double x, X) {
                os << x << " ";
            }
            os << std::endl;
        }
    }

    return os;
}
