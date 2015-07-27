#include "debruijn_graph.h"
#include "kmer_tbl.h"

#include <boost/foreach.hpp>

DeBruijnGraph::DeBruijnGraph(const KmerTable& tbl) {
    tbl.buildDeBruijn(this);
}

DeBruijnGraph::~DeBruijnGraph() {
}

void DeBruijnGraph::addKmer(const Kmer& kmer, size_t weight) {
    size_t k = kmer.length();
    Kmer key = kmer.subKmer(0, k - 1);
    Nucleotide::Code nucleotide = kmer.nucleotide(k - 1);

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
