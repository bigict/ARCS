#include "debruijn_graph.h"
#include "kmer_tbl.h"

#include <tr1/unordered_set>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.debruijn"));

DeBruijnGraph::DeBruijnGraph(const KmerTable& tbl) : _K(tbl.K()) {
    LOG4CXX_DEBUG(logger, boost::format("build debruijn graph from kmer hash table with K=[%d]") % _K);

    tbl.buildDeBruijn(this);
}

DeBruijnGraph::~DeBruijnGraph() {
}

void DeBruijnGraph::addKmer(const Kmer& kmer, size_t weight) {
    size_t k = kmer.length();
    Kmer key = kmer.subKmer(0, k - 1);
    Nucleotide::Code nucleotide = kmer.nucleotide(k - 1);

    LOG4CXX_TRACE(logger, boost::format("addKmer: %s with key=[%s], edge=[%c]") % kmer % key % Nucleotide::code2char(nucleotide));

    NodeList::iterator it = _nodelist.find(key);
    if (it != _nodelist.end()) {
        it->second.count[nucleotide] += weight;
    } else {
        _nodelist[key] = Node(nucleotide, 1);
    }
}

void DeBruijnGraph::removeKmer(const Kmer& kmer) {
    size_t k = kmer.length();
    Kmer key = kmer.subKmer(0, k - 1);
    Nucleotide::Code nucleotide = kmer.nucleotide(k - 1);

    LOG4CXX_TRACE(logger, boost::format("removeKmer: %s with key=[%s], edge=[%c]") % kmer % key % Nucleotide::code2char(nucleotide));
    
    NodeList::iterator it = _nodelist.find(key);
    if (it != _nodelist.end()) {
        if (it->second.count[nucleotide] > 0) {
            it->second.count[nucleotide] = 0;
        }
        if (!it->second) {
            _nodelist.erase(it);
        }
    }
}

void DeBruijnGraph::compact() {
    typedef std::tr1::unordered_map< Kmer, size_t, KmerHasher > CountingList;

    LOG4CXX_DEBUG(logger, boost::format("compact with %d nodes, begin") % _nodelist.size());

    CountingList merge_nodelist;
    {
    	LOG4CXX_DEBUG(logger, boost::format("counting indegree and outdegree begin"));
        // Counting in-degrees
        CountingList indegree, outdegree;

        for (NodeList::const_iterator it = _nodelist.begin(); it != _nodelist.end(); ++it) {
            for (size_t i = 0; i < _countof(it->second.count); ++i) {
                if (it->second.count[i] > 0) {
                    Kmer key = it->first.subKmer(1) + (Nucleotide::Code)i;
                    if (_nodelist.find(key) != _nodelist.end()) {
                        indegree[key]++;
                    }
                }
            }
            if (it->second.outdegree() == 1) {
                outdegree[it->first] = 1;
            }
        }

        for (CountingList::const_iterator it = outdegree.begin(); it != outdegree.end(); ++it) {
            if (indegree.find(it->first) != indegree.end()) {
                merge_nodelist[it->first] = 0;
            }
        }
    	LOG4CXX_DEBUG(logger, boost::format("counting indegree and outdegree end"));
    }

/*
    for (CountingList::iterator it = merge_nodelist.begin(); it != merge_nodelist.end(); ++it) {
    }
*/

    LOG4CXX_DEBUG(logger, boost::format("compact with %d nodes, end") % _nodelist.size());
}
