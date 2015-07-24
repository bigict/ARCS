#include "kmer_tbl.h"
#include "kseq.h"

#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.kmer_tbl"));

bool KmerTable::read(std::istream& stream) {
    if (!stream) {
        LOG4CXX_WARN(logger, boost::format("invalid input stream"));
        return false;
    }

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table begin"));

    DNASeqReader reader(stream);

    DNASeq seq;
    while (reader.read(seq)) {
    }

    LOG4CXX_DEBUG(logger, boost::format("construct kmer table end"));

    return true;
}
