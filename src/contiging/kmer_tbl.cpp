#include "kmer_tbl.h"
#include "debruijn_graph.h"
#include "kseq.h"

#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.kmer_tbl"));

KmerTable::KmerTable(size_t K, bool do_filter) : _K(K), _do_filter(do_filter) {
    BOOST_ASSERT(_K > 0);
}

KmerTable::~KmerTable() {
}

bool KmerTable::read(std::istream& stream) {
    if (!stream) {
        LOG4CXX_WARN(logger, boost::format("invalid input stream"));
        return false;
    }

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table begin"));

    DNASeqReader reader(stream);

    DNASeq read;
    while (reader.read(read)) {
        for (size_t i = 0,j = _K; j < read.seq.length(); ++i,++j) {
            Kmer kmer(read.seq, i, j);
            _hash_tbl[kmer]++;
        }
    }

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table end"));

    return true;
}

void KmerTable::buildDeBruijn(DeBruijnGraph* graph) const {
    LOG4CXX_DEBUG(logger, boost::format("construct de bruijn graph begin"));

    BOOST_ASSERT(graph != NULL);

    for (KmerList::const_iterator it = _hash_tbl.begin(); it != _hash_tbl.end(); ++it) {
        graph->addKmer(it->first, it->second);
    }

    LOG4CXX_DEBUG(logger, boost::format("construct de bruijn graph end"));
}
