#include "debruijn_graph.h"
#include "kmer_tbl.h"

#include <list>
#include <numeric>
#include <set>
#include <tr1/unordered_set>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.debruijn"));

DeBruijnGraph::DeBruijnGraph(const KmerTable& tbl) : _K(tbl.K()), _average(0) {
    LOG4CXX_DEBUG(logger, boost::format("build debruijn graph from kmer hash table with K=[%d]") % _K);

    tbl.statistics(&_average, NULL);
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
        if (it == _nodelist.end()) {
            // Add a virtual node
            Node n;
            n.parents[node] = 0;
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

bool DeBruijnGraph::findKmer(const Kmer& kmer, Node* node) const {
    return false;
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

void DeBruijnGraph::compact() {
    //typedef std::tr1::unordered_set< Kmer, KmerHasher > CountingList;
    typedef std::set< Kmer > CountingList;

    LOG4CXX_DEBUG(logger, boost::format("compact with %d nodes, begin") % _nodelist.size());

    // Xxxx
    CountingList merge_nodelist;
    {
        LOG4CXX_DEBUG(logger, boost::format("counting indegree and outdegree begin"));

        // Counting in-degrees & out-degrees
        for (NodeList::const_iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
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
        for (NodeList::const_iterator k = _nodelist.find(j->second.parents.begin()->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != j->first; k = _nodelist.find(k->second.parents.begin()->first)) {
            group.push_front(k->first);
        }

        // Look forward
        for (NodeList::const_iterator k = _nodelist.find(j->second.children.begin()->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != j->first; k = _nodelist.find(k->second.children.begin()->first)) {
            group.push_back(k->first);
        }

        // Merge nodes
        if (group.size() > 1) {
            Kmer key;

            LOG4CXX_DEBUG(logger, boost::format("compact key begin"));
            BOOST_FOREACH(const Kmer& kmer, group) {
                size_t overlap = key.overlap(kmer);
                BOOST_ASSERT(overlap != -1);
                key += kmer.subKmer(key.length() - overlap);
                LOG4CXX_DEBUG(logger, boost::format("compact: key=[%s]/[%s]") % key % kmer);
            }
            LOG4CXX_DEBUG(logger, boost::format("compact key end"));

            Node val;

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

struct KmerLengthPlus {
    size_t operator()(size_t l, const Kmer& r) const {
        return l + r.length();
    }
};

struct KmerCoveragePlus {
    size_t operator()(size_t l, const Kmer& r) const {
        return l + r.length();
    }
};

struct KmerRemover {
    KmerRemover(DeBruijnGraph::NodeList* nodelist) : _nodelist(nodelist) {
    } 

    void operator()(const Kmer& kmer) const {
        DeBruijnGraph::NodeList::iterator k = _nodelist->find(kmer);
        if (k != _nodelist->end()) {
            // Remove children
            for (DeBruijnGraph::EdgeList::const_iterator i = k->second.children.begin(); i != k->second.children.end(); ++i) {
                DeBruijnGraph::NodeList::iterator j = _nodelist->find(i->first);
                if (j != _nodelist->end()) {
                    j->second.parents.erase(kmer);
                }
            }
            k->second.children.clear();

            // Remove parents
            for (DeBruijnGraph::EdgeList::const_iterator i = k->second.parents.begin(); i != k->second.parents.end(); ++i) {
                DeBruijnGraph::NodeList::iterator j = _nodelist->find(i->first);
                if (j != _nodelist->end()) {
                    j->second.children.erase(kmer);
                }
            }
            k->second.parents.clear();
        }
    }

    typedef std::set< Kmer > NoiseList;
    NoiseList noiselist;

private:
    DeBruijnGraph::NodeList* _nodelist;
};


void DeBruijnGraph::removeNoise() {
    LOG4CXX_DEBUG(logger, boost::format("remove noise with %d nodes, begin") % _nodelist.size());

    KmerRemover remover(&_nodelist);

    for (NodeList::const_iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
        // Condensed node
        if (i->second.indegree() == 1 && i->second.outdegree() == 1) {
            continue;
        }

        for (EdgeList::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {

            std::list< Kmer > group;
            for (NodeList::const_iterator k = _nodelist.find(j->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != i->first; k = _nodelist.find(k->second.children.begin()->first)) {
                group.push_back(k->first);
            }

            size_t length = i->first.length() + std::accumulate(group.begin(), group.end(), 0, KmerLengthPlus());
            size_t coverage = j->second + std::accumulate(group.begin(), group.end(), 0, KmerCoveragePlus());
            if (i->second.indegree() == 0) {
                if (length < (group.size() + 1) * (_K - 1) + _K + 1) {
                    for_each(group.begin(), group.end(), remover);
                } else if (coverage < (group.size() + 1) * _average / 5.0) {
                    for_each(group.begin(), group.end(), remover);
                }
            } else if (length < (group.size() + 1) * (_K - 1) + _K + 1 && coverage < (group.size() + 1) * _average / 5.0) {
                for_each(group.begin(), group.end(), remover);
            } else if (coverage <= (group.size() + 1) * 3.0 || coverage < (group.size() + 1) * _average / 10.0) {
                for_each(group.begin(), group.end(), remover);
            }
        }
    }
    
    BOOST_FOREACH(const Kmer& kmer, remover.noiselist) {
        _nodelist.erase(kmer);
    }

    LOG4CXX_DEBUG(logger, boost::format("remove %d noise nodes, %d nodes left") % _nodelist.size() % remover.noiselist.size());
}

size_t DeBruijnGraph::Node::average() const {
    if (!children.empty()) {
        size_t n = 0;
        for (EdgeList::const_iterator it = children.begin(); it != children.end(); ++it) {
            n += it->second;
        }
        return n / children.size();
    }
    return 0;
}

