#include "component.h"

#include <boost/algorithm/string.hpp>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("remove_repeats.component"));

std::ostream& operator << (std::ostream& os, const Component& component) {
    os << boost::format("component %d") % component.id << std::endl;

    if (!component.items.empty()) {
        os << component.items[0].contig;
        for (size_t i = 1; i < component.items.size(); ++i) {
            os << boost::format(" %d") % component.items[i].contig;
        }
    }
    os << std::endl;
    if (component.items.size() > 1) {
        os << component.items[0].contig;
        for (size_t i = 1; i < component.items.size() - 1; ++i) {
            os << boost::format(" %d") % component.items[i].gap;
        }
    }
    os << std::endl;
    return os;
}

void Component::init(size_t id, const ContigList& contigs, const GapList& gaps) {
    BOOST_ASSERT(contigs.size() == gaps.size());

    this->id = id;

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

    static boost::regex reg(">component (\\d+)");

    if (_stream) {
        int state = eStart;
        std::string line;

        size_t id;
        Component::ContigList contigs;
        Component::GapList gaps;
        boost::char_separator<char> sep(" ,\t");

        while (std::getline(_stream, line)) {
            if (state == eStart) {
                boost::smatch what;
                //if (boost::algorithm::starts_with(line, ">")) {    
                if (boost::regex_match(line, what, reg)) {    
                    id = boost::lexical_cast< size_t >(what[1]);
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
                    component.init(id, contigs, gaps);
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
