#include "component.h"

#include <boost/algorithm/string.hpp>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("remove_repeats.component"));

void Component::init(const ContigList& contigs, const GapList& gaps) {
    BOOST_ASSERT(contigs.size() == gaps.size());

    items.clear();

    ContigList::const_iterator i = contigs.begin();
    GapList::const_iterator j = gaps.begin();
    for (; i != contigs.end() && j != gaps.end(); ++i,++j) {
        items.push_back(Item(*i, *j));
    }
}

bool ComponentReader::read(Component& component) {
    enum {
        eStart,
        eId,
        eGap,
    };

    if (_stream) {
        int state = eStart;
        std::string line;

        Component::ContigList contigs;
        Component::GapList gaps;
        boost::char_separator<char> sep(" ,\t");

        while (std::getline(_stream, line)) {
            if (state == eStart) {
                if (boost::algorithm::starts_with(line, ">")) {    
                    state = eId;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fa=>invalid line for not start with >: %s") % line);
                    return false;
                }
            } else if (state == eId) {
                boost::tokenizer< boost::char_separator< char > > toker(line, sep);
                BOOST_FOREACH(const std::string& id, toker)  {
                    contigs.push_back(boost::lexical_cast< size_t >(id));
                }
                state = eGap;
            } else if (state == eGap) {
                boost::tokenizer< boost::char_separator< char > > toker(line, sep);
                BOOST_FOREACH(const std::string& gap, toker)  {
                    gaps.push_back(boost::lexical_cast< int >(gap));
                }
                gaps.push_back(0);//finall contige 

                if (contigs.size() == gaps.size()) {
                    component.init(contigs, gaps);
                    return true;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fa=>invalid line for id size != gap size: %s") % line);
                    return false;
                }
            }
        }
    }

    return false;
}
