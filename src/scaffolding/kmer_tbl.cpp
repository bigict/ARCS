#include "kmer_tbl.h"

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffolding.kmer_tbl"));

void KmerTable::addKmer( const Kmer& o, std::pair<size_t, long> pos) {
    _kmerlist.insert(std::make_pair( o, pos ));
}

std::pair<size_t, long> KmerTable::findPos(const Kmer& o) const {
    KmerList::const_iterator it = _kmerlist.find(o);
    if (it != _kmerlist.end()) {
        return it->second;
    }
    return std::make_pair(-1, -1);
}

std::ostream& operator<<(std::ostream& os, const KmerTable& tbl) {
	os << boost::format("kmer list size: %d") % tbl._kmerlist.size() << std::endl;

	for (KmerTable::KmerList::const_iterator it = tbl._kmerlist.begin(); it != tbl._kmerlist.end(); ++it) {
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

size_t BuildKmerTable(size_t K, size_t insert_size, const ContigList& contigs, const ComponentList& components, KmerList& tbl) {
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
