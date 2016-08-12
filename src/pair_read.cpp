#include "pair_read.h"
#include "kmer.h"
#include "kseq.h"

#include <fstream>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.PairRead"));

class PairReadReader {
public:
    PairReadReader(std::istream& stream1, std::istream& stream2) : _reader1(stream1), _reader2(stream2) {
    }

    bool read(PairRead& read) {
        DNASeq read1, read2;
        if (_reader1.read(read1) && _reader2.read(read2)) {
            read.set1 = read1.seq;
            read.set2 = read2.seq;
            return true;
        }
        return false;
    }
private:
    DNASeqReader _reader1;
    DNASeqReader _reader2;
};

bool ReadPairReads(std::istream& stream1, std::istream& stream2, PairReadList& pair_reads) {
    PairReadReader reader(stream1, stream2);
    PairRead pair_read;
    while (reader.read(pair_read)) {
        pair_reads.push_back(pair_read);
    }
    return true;
}

bool ReadSinglePairReadsFile(const std::string& file1, const std::string& file2, PairReadList& pair_reads) {
    std::ifstream stream1(file1.c_str()), stream2(file2.c_str());
    return ReadPairReads(stream1, stream2, pair_reads);
}

bool ReadPairReads(const std::string& files1, const std::string& files2, PairReadList& pair_reads) {
    std::vector<std::string> r1;
    boost::split(r1, files1, boost::is_any_of(";:"));
    std::vector<std::string> r2;
    boost::split(r2, files2, boost::is_any_of(";:"));
    BOOST_ASSERT(r1.size() != 0 && r1.size() == r2.size());
    for(size_t i = 0; i < r1.size(); ++i) {
        LOG4CXX_DEBUG(logger, boost::format("read file1: %s read file2: %s") % r1[i] % r2[i]);
        if(!ReadSinglePairReadsFile(r1[i] , r2[i], pair_reads)) {
            return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const PairReadList& pair_reads) {
    BOOST_FOREACH(const PairRead& pair_read, pair_reads) {
        os << boost::format("%s\t%s") % pair_read.set1 % pair_read.set2 << std::endl; 
    }
    return os;
}

