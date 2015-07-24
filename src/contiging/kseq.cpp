#include "kseq.h"

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.kseq"));

int DNASeq::quality() const {
    return 0;
}

bool DNASeqReader::read(DNASeq& sequence) {
    enum {
        kName, 
        kSequence, 
        kName2, 
        kQuality, 
    };

    if (_stream) {
        int state = kName;
        std::string buf;

        while (std::getline(_stream, buf)) {
            if (state == kName) {
                if (boost::algorithm::starts_with(buf, "@")) {
                    sequence.name = buf.substr(1);
                    state = kSequence;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fastq=>invalid line for sequence name: %s") % buf);
                    return false;
                }
            } else if (state == kSequence) {
                sequence.seq = buf;
                state = kName2;
            } else if (state == kName2) {
                if (boost::algorithm::starts_with(buf, "+") && boost::algorithm::ends_with(buf, sequence.name)) {
                    state = kQuality;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fastq=>names aren't equal: %s") % buf);
                    return false;
                }
            } else if (state == kQuality) {
                if (buf.length() == sequence.seq.length()) {
                    sequence._quality = buf;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fastq=>length of sequence and quality are not equal: %s") % buf);
                    return false;
                }
            }
        }
    }

    return false;
}
