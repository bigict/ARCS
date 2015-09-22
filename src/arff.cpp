#include "arff.h"

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.ARFFReader"));

bool ARFFReader::read(ARFF& arff) {
    LOG4CXX_DEBUG(logger, "read ARFF from stream begin");

    if (!_stream) {
        LOG4CXX_WARN(logger, "read ARFF from stream failed");
        return false;
    }

    enum {
        kAttribute = 0, 
        kData, 
        kTuple, 
    };
    int state = kAttribute;

    static boost::regex mapping[3] = {
        boost::regex("@attribute\\s+([^\\s]+)\\s+(.+)"), 
        boost::regex("@data\\s*"), 
        boost::regex("([ACGTN]+)\\t(\\d+)")
    };
    static boost::char_separator< char > sep("\t");
    boost::smatch what;

    std::string line;
    while (std::getline(_stream, line)) {
        boost::algorithm::trim(line);
        if (line.empty()) continue;

        if (state == kAttribute) {
            if (boost::regex_match(line, what, mapping[kAttribute])) {
                arff.attrs[what[1]] = what[2];
            } else if (boost::regex_match(line, what, mapping[kData])) {
                state = kTuple;
            } else {
                return false;
                LOG4CXX_WARN(logger, "read ARFF from stream failed@attribute");
            }
        } else if (state == kTuple) {
            ARFF::ItemPtr item = new ARFF::Item(_num_fields);
            boost::tokenizer< boost::char_separator< char > > tok(line, sep);
            size_t k = 0;
            BOOST_FOREACH(const std::string& i, tok) {
                (*item)[k++] = i;
            }
            arff.data.push_back(item);
        }
    }
    
    LOG4CXX_DEBUG(logger, "read ARFF from stream end");

    return true;
}
