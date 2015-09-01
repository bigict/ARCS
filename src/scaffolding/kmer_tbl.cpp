#include "kmer_tbl.h"

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffolding.kmer_tbl"));

std::ostream& operator<<(std::ostream& os, const KmerList& tbl) {
	os << boost::format("kmer list size: %d") % tbl.size() << std::endl;

	for (KmerList::const_iterator it = tbl.begin(); it != tbl.end(); ++it) {
		os << boost::format("%d\t%s\t%d\t%d") % it->first.length() % it->first.sequence() % it->second.first % it->second.second << std::endl;
    }

	return os;
}

size_t _BuildKmerTable(size_t K, const ContigList& contigs, size_t component_no, const Component& component, KmerList& tbl) {
    size_t num = 0;
    
    long idx = 0;
    for (size_t k = 0; k < component._contig_id.size(); ++k) {
        const Contig& contig = contigs[component._contig_id[k]];
        for (size_t i = 0,j = K; j <= contig.seq.length(); ++i,++j) {
            Kmer kmer(contig.seq, i, j);
            tbl.insert(std::make_pair(kmer, KmerPosition(component_no, idx)));
            ++idx;
            ++num;
        }
        idx += component._gap[k];
    }

    return num;
}

size_t _BuildKmerTable(size_t K, size_t insert_size, const ContigList& contigs, size_t component_no, const Component& component, KmerList& tbl) {
    size_t num = 0;
    
    size_t cutoff = 2 * insert_size;
    long idx = 0;
    if (component._contig_id.size() == 1) {
        const  Contig& contig = contigs[component._contig_id[0]];
        for (size_t i = 0, j = K; j <= contig.seq.length(); ++i,++j) {
            if (idx <= cutoff || idx >= component.length() - 1 - K + 1 - cutoff) {
                Kmer kmer(contig.seq, i, j);
                tbl.insert(std::make_pair(kmer, KmerPosition(component_no, idx)));
                ++num;
            }
            ++idx;
        }
    } else {
        for (size_t k = 0; k < component._contig_id.size(); ++k) {
            const Contig& contig = contigs[component._contig_id[k]];
            for (size_t i = 0,j = K; j <= contig.seq.length(); ++i,++j) {
                if (idx <= cutoff || idx >= component.length() - 1 - cutoff) {
                    Kmer kmer(contig.seq, i, j);
                    tbl.insert(std::make_pair(kmer, KmerPosition(component_no, idx)));
                    ++num;
                }
                ++idx;
            }
            idx += component._gap[k];
        }
    }

    return num;
}

size_t BuildKmerTable_insertsize(size_t K, size_t insert_size, const ContigList& contigs, const ComponentList& components, KmerList& tbl) {
    size_t idx = 0, num = 0;

    LOG4CXX_DEBUG(logger, boost::format("build kmer table for insert size begin"));
    BOOST_FOREACH(const Component& component, components) {
        if (component._contig_id.empty() || component.length() < 2 * insert_size) continue;
        num += _BuildKmerTable(K, contigs, idx, component, tbl);
        ++idx;
    }
    LOG4CXX_DEBUG(logger, boost::format("build kmer table for insert size end"));

    return num;
}

size_t BuildKmerTable_pairends(size_t K, size_t insert_size, size_t edge_cutoff, const ContigList& contigs, const ComponentList& components, KmerList& tbl) {
    size_t idx = 0, num = 0;

    LOG4CXX_DEBUG(logger, boost::format("build kmer table for pair ends begin"));
    BOOST_FOREACH(const Component& component, components) {
        if (component._contig_id.empty() || (component._contig_id.size() == 1 && contigs[component._contig_id[0]].seq.length() < edge_cutoff)) {
            ++idx;
            continue;
        }
        num += _BuildKmerTable(K, insert_size, contigs, idx, component, tbl);
        ++idx;
    }
    LOG4CXX_DEBUG(logger, boost::format("build kmer table for pair ends end"));

    return num;
}

