#include "kseq.h"

#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.kseq"));

void DNASeq::make_complement() {
	std::reverse(seq.begin(), seq.end());
	std::reverse(quality.begin(), quality.end());

	static std::map< char, char > mapping = boost::assign::map_list_of
		('A', 'T')
		('C', 'G')
		('G', 'C')
		('T', 'A')
		('N', 'N');

	for (size_t i = 0; i < seq.size(); ++i) {
		seq[i] = mapping[seq[i]];
	}
}

std::ostream& operator << (std::ostream& os, const DNASeq& seq) {
    os << '@' << seq.name << std::endl;
    os << seq.seq << std::endl;
    os << '+' << std::endl;
    os << seq.quality << std::endl;
    return os;
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
                if (boost::algorithm::starts_with(buf, "+") && (buf.length() == 1 || boost::algorithm::ends_with(buf, sequence.name))) {
                    state = kQuality;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fastq=>names aren't equal: %s") % buf);
                    return false;
                }
            } else if (state == kQuality) {
                if (buf.length() == sequence.seq.length()) {
                    sequence.quality = buf;
                    cutoff(sequence);
                    return true;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fastq=>length of sequence and quality are not equal: %s") % buf);
                    return false;
                }
            }
        }
    }

    return false;
}

void DNASeqReader::cutoff(DNASeq& sequence) const {
    size_t length = sequence.seq.length();
    if (_percent < 1.0) {
        length = (size_t)(length * _percent);
    }
    length = std::min(length, _read_cutoff);
    if (sequence.seq.length() > length) {
        sequence.seq = sequence.seq.substr(0, length);
        sequence.quality = sequence.quality.substr(0, length);
    }
}
