#include "debruijn_graph.h"
#include "kmer_tbl.h"

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

    tbl.statistics(&_average, NULL);

    LOG4CXX_DEBUG(logger, boost::format("average coverage=[%f]") % _average);

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
        return l + r.length();
    }

    static size_t length(size_t K, const Kmer& kmer) {
        return kmer.length() - (K-1) + 1;
    }
private:
    size_t _K;
};

struct KmerCoveragePlus {
    KmerCoveragePlus(DeBruijnGraph::NodeList* nodelist, size_t K) : _nodelist(nodelist), _K(K) {
    } 

    size_t operator()(size_t l, const Kmer& r) const {
        DeBruijnGraph::NodeList::const_iterator i = _nodelist->find(r);
        if (i != _nodelist->end()) {
            for (DeBruijnGraph::EdgeList::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
                l += coverage(_K, j->first, j->second);
            }
        }
        //return l + _nodelist->find(r)->second.parents.begin()->second;
    }

    typedef std::set< Kmer > NoiseList;

    static size_t coverage(size_t K, const Kmer& kmer, size_t n) {
        return n * KmerLengthPlus::length(K, kmer);
    }
private:
    DeBruijnGraph::NodeList* _nodelist;
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
        for (NodeList::const_iterator k = _nodelist.find(j->second.parents.begin()->first); k != _nodelist.end() && k->second.indegree() <= 1 && k->second.outdegree() == 1 && k->first != j->first; k = (k->second.parents.empty() ? _nodelist.end() : _nodelist.find(k->second.parents.begin()->first))) {
            group.push_front(k->first);
        }

        // Look forward
        for (NodeList::const_iterator k = _nodelist.find(j->second.children.begin()->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != j->first; k = (k->second.children.empty() ? _nodelist.end() : _nodelist.find(k->second.children.begin()->first))) {
            group.push_back(k->first);
        }

        // Merge nodes
        if (group.size() > 1) {
            Kmer key;

            BOOST_FOREACH(const Kmer& kmer, group) {
                size_t overlap = key.overlap(kmer);
                BOOST_ASSERT(overlap != -1);
                key += kmer.subKmer(key.length() - overlap);//key->kmer????
                LOG4CXX_TRACE(logger, boost::format("compact: key=[%s]/[%s]") % key % kmer);
            }

            Node val;

            // Average coverage
            size_t length = std::accumulate(group.begin(), group.end(), 0, KmerLengthPlus(_K));
            size_t coverage = std::accumulate(group.begin(), group.end(), 0, KmerCoveragePlus(&_nodelist, _K));

            // Set parrents
            Node& front(_nodelist[group.front()]);
            val.parents = front.parents;
            for (EdgeList::const_iterator it = front.parents.begin(); it != front.parents.end(); ++it) {
                NodeList::iterator k = _nodelist.find(it->first);
                if (k != _nodelist.end()) {
                    //k->second.children[key] = coverage / group.size(); // Use average coverage here.
                    k->second.children[key] = k->second.children[group.front()];
                    k->second.children.erase(group.front());
                }
            }

            // Set children
            Node& back(_nodelist[group.back()]);
            val.children = back.children;
            for (EdgeList::iterator it = back.children.begin(); it != back.children.end(); ++it) {

                it->second = (coverage + length / 2) / length; // Use average coverage here.

                NodeList::iterator k = _nodelist.find(it->first);
                if (k != _nodelist.end()) {
                    //k->second.parents[key] = coverage / group.size();
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
#define KMER_REMOVER_TRANNSACTION_BEGIN(remover)     remover.transaction_begin();
#define KMER_REMOVER_TRANNSACTION_END(remover)       remover.transaction_end();
struct KmerRemover {
    KmerRemover(DeBruijnGraph::NodeList* nodelist) : _nodelist(nodelist) {
    } 

    void operator()(const Kmer& kmer) {
        _translist.insert(kmer);
    }

    void transaction_begin() {
        _translist.clear();
    }
    void transaction_end() {
        for_each(_translist.begin(), _translist.end(), boost::bind(&KmerRemover::clear, this, _1));
    }

    void clear(const Kmer& kmer) {
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
            noiselist.insert(kmer);
        }
    }

    typedef std::set< Kmer > NoiseList;
    NoiseList noiselist;

private:
    NoiseList _translist;
    DeBruijnGraph::NodeList* _nodelist;
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

        KMER_REMOVER_TRANNSACTION_BEGIN(remover);

        for (EdgeList::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {

            std::list< Kmer > group;
            // Look forward
            for (NodeList::const_iterator k = _nodelist.find(j->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != i->first; k = _nodelist.find(k->second.children.begin()->first)) {
                group.push_back(k->first);
            }

            // front && back
            NodeList::const_iterator front(i), back(group.empty() ? _nodelist.find(j->first) : _nodelist.find(_nodelist.find(group.back())->second.children.begin()->first));

            // length && coverage
            size_t length = KmerLengthPlus::length(_K, front->first) + std::accumulate(group.begin(), group.end(), 0, KmerLengthPlus(_K)) + KmerLengthPlus::length(_K, back->first);
            size_t coverage = (KmerCoveragePlus::coverage(_K, j->first, j->second) + std::accumulate(group.begin(), group.end(), 0, KmerCoveragePlus(&_nodelist, _K)) + length / 2) / length; 
            LOG4CXX_DEBUG(logger, boost::format("length=[%d],avg_coverage=[%d]") % length % coverage);

            // rules
            if (front->second.indegree() == 0 || back->second.outdegree() == 0 ) { // tips
                if (_K - 1 + length - 1 < 2 * _K || coverage < _average / 5.0) {
                    if (front->second.indegree() == 0) {
                        remover(front->first);
                    } else {
                        remover(back->first);
                    }
                    for_each(group.begin(), group.end(), remover);
                }
            } else if ((_K - 1 + length - 1 < 2 * _K && coverage < _average / 5.0) || (coverage < _average / 10.0)) { // buble or link ?????????
                if (front->second.indegree() == 0) {
                    remover(front->first);
                } else {
                    remover(back->first);
                }
                for_each(group.begin(), group.end(), remover);
            }
        }

        KMER_REMOVER_TRANNSACTION_END(remover);
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

std::ostream& operator << (std::ostream& os, const DeBruijnGraph& graph) {
    DeBruijnGraph::EdgeList index;

    size_t num = 1;
    for (DeBruijnGraph::NodeList::const_iterator i=graph._nodelist.begin(); i!=graph._nodelist.end(); ++i) {
        if ( i->second.indegree() != 0 || i->second.outdegree() != 0){
            index[ i->first ] = ++num;
        }
    }
    for (DeBruijnGraph::NodeList::const_iterator i=graph._nodelist.begin(); i!=graph._nodelist.end(); ++i) {
        for(DeBruijnGraph::EdgeList::const_iterator j=i->second.children.begin(); j!=i->second.children.end(); ++j) {
            Kmer key = i->first;
            size_t overlap = key.overlap( j->first );
            key += j->first.subKmer(key.length() - overlap);
            //os << overlap << " " << j->first.subKmer(j->first.length() - overlap) << std::endl;
            os << index[ i->first ] << "\t" << index[ j->first] << "\t" << key << std::endl;
            os << j->second << std::endl;
        }
    }
    return os;
}

void DeBruijnGraph::outputCdbgGraphAndComponent0(const std::string& cdbg_file, const std::string& component0_file) const {
    
    std::ofstream cdbg(cdbg_file.c_str());
    std::ofstream comp(component0_file.c_str());
    if (!cdbg || !comp) {
        LOG4CXX_DEBUG(logger, boost::format("open write file error"));
        exit(1);
    }
    size_t contig_num = 0;
    size_t unique_num = 0;
    for (NodeList::const_iterator i=_nodelist.begin(); i!=_nodelist.end(); ++i) {
        for(EdgeList::const_iterator j=i->second.children.begin(); j!=i->second.children.end(); ++j) {
            size_t copy_num = size_t((j->second + _average/2) / _average);
            Kmer key = i->first;
            size_t overlap = key.overlap( j->first );
            key += j->first.subKmer(key.length() - overlap);
            //os << overlap << " " << j->first.subKmer(j->first.length() - overlap) << std::endl;
            cdbg << ">seq_" << contig_num << "\t" << copy_num << std::endl;
            cdbg << key << std::endl;
            if (copy_num <= 1) {
                comp << ">component " << "\t" << unique_num << std::endl;
                comp << contig_num << std::endl << std::endl;
                ++ unique_num;
            }
            ++ contig_num;
        }
    }
}


