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
    size_t operator()(size_t l, const Kmer& r) const {
        return l + r.length();
    }
};

struct KmerCoveragePlus {
    KmerCoveragePlus(DeBruijnGraph::NodeList* nodelist) : _nodelist(nodelist) {
    } 

    size_t operator()(size_t l, const Kmer& r) const {
        return l + _nodelist->find(r)->second.parents.begin()->second;
    }

    typedef std::set< Kmer > NoiseList;

private:
    DeBruijnGraph::NodeList* _nodelist;
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

            BOOST_FOREACH(const Kmer& kmer, group) {
                size_t overlap = key.overlap(kmer);
                BOOST_ASSERT(overlap != -1);
                key += kmer.subKmer(key.length() - overlap);//key->kmer????
                LOG4CXX_TRACE(logger, boost::format("compact: key=[%s]/[%s]") % key % kmer);
            }

            Node val;

            // Average coverage
            size_t coverage = std::accumulate(group.begin(), group.end(), 0, KmerCoveragePlus(&_nodelist));

            // Set parrents
            Node& front(_nodelist[group.front()]);
            val.parents = front.parents;
            for (EdgeList::const_iterator it = front.parents.begin(); it != front.parents.end(); ++it) {
                NodeList::iterator k = _nodelist.find(it->first);
                if (k != _nodelist.end()) {
                    k->second.children[key] = coverage / group.size(); // Use average coverage here.
                    k->second.children.erase(group.front());
                }
            }

            // Set children
            Node& back(_nodelist[group.back()]);
            val.children = back.children;
            for (EdgeList::iterator it = back.children.begin(); it != back.children.end(); ++it) {
                NodeList::iterator k = _nodelist.find(it->first);
                if (k != _nodelist.end()) {
                    k->second.parents[key] = coverage / group.size();
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

    void operator()(const Kmer& kmer) const {
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
        }
    }

    typedef std::set< Kmer > NoiseList;
    NoiseList noiselist;

private:
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

        for (EdgeList::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {

            std::list< Kmer > group;
            // Look forward
            for (NodeList::const_iterator k = _nodelist.find(j->first); k != _nodelist.end() && k->second.indegree() == 1 && k->second.outdegree() == 1 && k->first != i->first; k = _nodelist.find(k->second.children.begin()->first)) {
                group.push_back(k->first);
            }

            size_t length = i->first.length() + std::accumulate(group.begin(), group.end(), 0, KmerLengthPlus()) + (group.size()!=0 ? _nodelist.find(*--group.end())->second.children.begin()->first.length() : j->first.length());
            //size_t coverage = j->second + std::accumulate(group.begin(), group.end(), 0, KmerCoveragePlus());
            size_t coverage = j->second * (j->first.length() - _K + 2);
            size_t cov_num = j->first.length() -_K + 2;
            for (std::list< Kmer >::const_iterator k=group.begin(); k!=group.end(); ++k) {
                NodeList::const_iterator it = _nodelist.find( *k );
                BOOST_ASSERT(it != _nodelist.end());
                size_t kmer_len = it->second.children.begin() -> first.length() - _K + 2;
                coverage += it->second.children.begin() -> second * kmer_len;
                            
                cov_num += kmer_len;
            }
            size_t avg_coverage = (coverage + cov_num / 2 ) / cov_num;

            if (group.size() != 0) {
                std::cout << *group.begin() << std::endl;
            }
            std::cout << "length:" << length << " avg_cov:" << avg_coverage << std::endl;
            if (group.size() == 0) {
                //group.push_back( j->first ); 
                if (i->second.indegree() == 0 || _nodelist.find(j->first)->second.outdegree() ) {
                    if (length - (group.size() + 1) * (_K - 2) < 2 * _K) {
                        _nodelist.find(i->first)->second.children.erase( j->first );
                        _nodelist.find(j->first)->second.parents.erase( i->first );
                        //for_each(group.begin(), group.end(), remover);
                    } else if (avg_coverage < _average / 5.0) {
                        _nodelist.find(i->first)->second.children.erase( j->first );
                        _nodelist.find(j->first)->second.parents.erase( i->first );
                        //for_each(group.begin(), group.end(), remover);
                    }
                } else if (length - (group.size() + 1) * (_K - 2) < 2 * _K && avg_coverage < _average / 5.0) {
                    _nodelist.find(i->first)->second.children.erase( j->first );
                    _nodelist.find(j->first)->second.parents.erase( i->first );
                    //for_each(group.begin(), group.end(), remover);
                } else if (length <= (group.size() + 1) * 3.0 || avg_coverage < _average / 10.0) {
                    _nodelist.find(i->first)->second.children.erase( j->first );
                    _nodelist.find(j->first)->second.parents.erase( i->first );
                    //for_each(group.begin(), group.end(), remover);
                }
                
            }
            else {
                if (i->second.indegree() == 0 || _nodelist.find( _nodelist.find(*--group.end())->second.children.begin()->first)->second.outdegree() == 0) {
                    if (length - (group.size() + 1) * (_K - 2) < 2 * _K) {
                        for_each(group.begin(), group.end(), remover);
                    } else if (avg_coverage < _average / 5.0) {
                        for_each(group.begin(), group.end(), remover);
                    }
                } else if (length - (group.size() + 1) * (_K - 2) < 2 * _K && avg_coverage < _average / 5.0) {
                    for_each(group.begin(), group.end(), remover);
                } else if (length <= (group.size() + 1) * 3.0 || avg_coverage < _average / 10.0) {
                    for_each(group.begin(), group.end(), remover);
                }
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
            os << index[ i->first ] << "\t" << index[ j->first] << "\t" << key << std::endl;
            os << j->second << std::endl;
        }
    }
    return os;
}


