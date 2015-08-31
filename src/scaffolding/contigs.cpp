#include "contigs.h"

#include <fstream>
#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>


#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffolding.contig"));

std::ostream& operator<<(std::ostream& os, const ContigList& contigs) {
    size_t index = 0;
    BOOST_FOREACH(const Contig& contig, contigs) {
        os << index << "\t" << contig.seq << "\t" << contig.copy_num << "\t" << contig.seq.length() << std::endl;
        ++ index; 
    }
    return os;
} 

bool ContigReader::read(Contig& contig) {
    enum {
        eCopy,
        eContig, 
    };

    if (_stream) {
        static boost::regex reg(">seq_(\\d+)[ \t]+(\\d+)");

        int state = eCopy;
        std::string buf;

        while (std::getline(_stream, buf)) {
            if (state == eCopy) {
                boost::smatch what;
                if (boost::regex_match(buf, what, reg)) {
                    //contig.id = boost::lexical_cast< size_t >(what[1]);
                    contig.copy_num = boost::lexical_cast< int > (what[2]);
                    state = eContig;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fa=>invalid line for edge name:not start with > %s") % buf);
                    return false;
                }
            } else if (state == eContig) {
                contig.seq = buf;
                state = eCopy;
                return true;
            }
        }
    }
    return false;
}

bool ReadContigs(std::istream& stream, ContigList& contigs) {
    if (!stream) {
        LOG4CXX_ERROR(logger, boost::format("can not read contig from stream"));
        return false;
    }
    LOG4CXX_DEBUG(logger, boost::format("read contig from stream begin"));

    ContigReader reader(stream);
    Contig contig;
    while (reader.read(contig)) {
        contigs.push_back(contig);
    }

    LOG4CXX_DEBUG(logger, boost::format("read contig from stream end"));

    return true;
}

bool ReadContigs(const std::string& file, ContigList& contigs) {
    std::ifstream stream(file.c_str());
    return ReadContigs(stream, contigs);
}

size_t GenomeLength(const ContigList& contigs, size_t K) {
    size_t length = 0;
    BOOST_FOREACH(const Contig& contig, contigs) {
        length += contig.seq.size() -K + 1;
    }
    return length;
}

