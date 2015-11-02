#include "condensed_debruijn_graph_reader.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.CondensedDeBruijnGraph"));

bool CondensedDeBruijnGraphReader::read(CondensedDeBruijnGraphEdge& edge) {
    enum {
        kEdge, 
        kCoverage
    };
    if (_stream) {
        static boost::regex reg("(\\d+)\\s+(\\d+)\\s+([ACGTN]+)[\\s+.*]*");

        try {
            int state = kEdge;
            std::string line;
            while (std::getline(_stream, line)) {
                boost::algorithm::trim(line);
                if (line.empty()) continue;
                if (state == kEdge) {
                    boost::smatch what;
                    if (boost::regex_match(line, what, reg)) {
                        edge.i = boost::lexical_cast< size_t >(what[1]);
                        edge.j = boost::lexical_cast< size_t >(what[2]);
                        edge.seq = what[3];
                        state = kCoverage;
                    } else {
                        return false;
                    }
                } else if (state == kCoverage) {
                    edge.coverage = boost::lexical_cast< double >(line);
                    return true;
                }
            }
        } catch (...) {
        }
    }
    return false;
}

