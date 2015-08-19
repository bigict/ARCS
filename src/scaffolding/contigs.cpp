#include "contigs.h"

#include <fstream>
#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scanfording.contig"));

ContigSet::ContigSet(std::istream& stream, size_t K) : _K(K) {
    init(stream);
}

ContigSet::ContigSet(const std::string& filename, size_t K) : _K(K) {
    std::ifstream stream(filename.c_str());
    init(stream);
}

void ContigSet::init(std::istream& stream) {
    LOG4CXX_INFO(logger, boost::format("read contig file begin"));

    if (!stream) {
        LOG4CXX_ERROR(logger,boost::format("can not find contig file"));
    }

    GENOME_LEN = 0;

    ContigReader reader(stream);
    Contig contig;
    while (reader.read(contig)) {
        GENOME_LEN += contig.contig.size() - _K + 1;
        contigs.push_back(contig);  
    }

    LOG4CXX_INFO(logger, boost::format("all contig length GENOME_LEN: %d") % GENOME_LEN);

    LOG4CXX_INFO(logger, boost::format("read contig file end"));
}

std::ostream& operator<<(std::ostream& os, const ContigSet& contigs) {
    size_t index = 0;
    BOOST_FOREACH(const Contig& contig, contigs.contigs) {
        os << index << "\t" << contig.contig << "\t" << contig.copy_num << "\t" << contig.contig.length() << std::endl;
        //os << boost::format(">seq_%d \t%d") % contig.id % contig.copy_num << std::endl;
        //os << contig.contig << std::endl;
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
                    contig.id = boost::lexical_cast< size_t >(what[1]);
                    contig.copy_num = boost::lexical_cast< int > (what[2]);
                    state = eContig;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fa=>invalid line for edge name:not start with > %s") % buf);
                    return false;
                }
            } else if (state == eContig) {
                contig.contig = buf;
                state = eCopy;
                return true;
            }
        }
    }
    return false;
}
