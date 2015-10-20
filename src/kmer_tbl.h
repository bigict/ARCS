#ifndef kmer_tbl_h_
#define kmer_tbl_h_

#include "component.h"
#include "contigs.h"
#include "kmer.h"

#include <iostream>
#include <tr1/unordered_map>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

typedef std::pair< size_t, long > KmerPosition;

template< size_t K >
size_t _BuildKmerTable_(size_t L, const ContigList& contigs, size_t component_no, const Component& component, KmerTable< K, KmerPosition >& tbl) {
    size_t num = 0;
    
    long idx = 0;
    for (size_t k = 0; k < component.contigs.size(); ++k) {
        const Contig& contig = contigs[component.contigs[k]];
        for (size_t i = 0,j = L; j <= contig.seq.length(); ++i,++j) {
            Kmer< K > kmer(contig.seq, i, j);
            tbl.insert(std::make_pair(kmer, KmerPosition(component_no, idx)));
            ++idx;
            ++num;
        }
        idx += component.gaps[k];
    }

    return num;
}

template< size_t K >
size_t _BuildKmerTable_(size_t L, size_t insert_size, const ContigList& contigs, size_t component_no, const Component& component, KmerMultiTable< K, KmerPosition >& tbl) {
    size_t num = 0;
    
    size_t cutoff = 2 * insert_size;
    long idx = 0;
    if (component.contigs.size() == 1) {
        const  Contig& contig = contigs[component.contigs[0]];
        for (size_t i = 0, j = L; j <= contig.seq.length(); ++i,++j) {
            if (idx <= cutoff || (component.length() >= L + cutoff && idx >= component.length() - L + 1 - 1 - cutoff)) {
                Kmer< K > kmer(contig.seq, i, j);
                tbl.insert(std::make_pair(kmer, KmerPosition(component_no, idx)));
                ++num;
            }
            ++idx;
        }
    } else {
        for (size_t k = 0; k < component.contigs.size(); ++k) {
            const Contig& contig = contigs[component.contigs[k]];
            for (size_t i = 0,j = L; j <= contig.seq.length(); ++i,++j) {
                if (idx <= cutoff || (component.length() >= L + cutoff && idx >= component.length() - L + 1 - 1 - cutoff)) {
                    Kmer< K > kmer(contig.seq, i, j);
                    tbl.insert(std::make_pair(kmer, KmerPosition(component_no, idx)));
                    ++num;
                }
                ++idx;
            }
            idx += component.gaps[k];
        }
    }

    return num;
}

template< size_t K >
size_t BuildKmerTable_insertsize(size_t L, size_t insert_size, const ContigList& contigs, const ComponentList& components, KmerTable< K, KmerPosition >& tbl) {
    size_t idx = 0, num = 0;

    BOOST_FOREACH(const Component& component, components) {
        if (component.contigs.empty() || component.length() < 2 * insert_size) continue;
        num += _BuildKmerTable_< K >(L, contigs, idx, component, tbl);
        ++idx;
    }

    return num;
}

template< size_t K >
size_t BuildKmerTable_pairends(size_t L, size_t insert_size, size_t edge_cutoff, const ContigList& contigs, const ComponentList& components, KmerMultiTable< K, KmerPosition >& tbl) {
    size_t idx = 0, num = 0;

    BOOST_FOREACH(const Component& component, components) {
        if (component.contigs.empty() || (component.contigs.size() == 1 && contigs[component.contigs[0]].seq.length() < edge_cutoff)) {
            ++idx;
            continue;
        }
        num += _BuildKmerTable_< K >(L, insert_size, contigs, idx, component, tbl);
        ++idx;
    }

    return num;
}

template< size_t K >
std::ostream& operator<<(std::ostream& os, const KmerTable< K, KmerPosition >& tbl) {
	os << boost::format("kmer list size: %d") % tbl.size() << std::endl;

	for (typename KmerTable< K, KmerPosition >::const_iterator it = tbl.begin(); it != tbl.end(); ++it) {
		os << boost::format("%d\t%s\t%d\t%d") % it->first.length() % it->first.sequence() % it->second.first % it->second.second << std::endl;
    }

	return os;
}

#endif // kmer_tbl_h_

