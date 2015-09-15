#include "component.h"

#include <fstream>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffolding.component"));

void Component::length(size_t K, const ContigList& contig_list) {
    _length = 0;
    if (!contigs.empty()) {
        _length = contig_list[contigs[0]].seq.length() + 1;//can change
        for (size_t i = 1; i < contig_list.size(); ++i) {
            _length += gaps[i - 1];
            BOOST_ASSERT(contig_list[contigs[i]].seq.size() >= K - 1);
            _length += contig_list[contigs[i]].seq.size() - K + 1;
        }
    }
}

bool ComponentReader::read(Component& component) {
    enum {
        eStart,
        eId,
        eGap,
    };

    static boost::regex reg(">component[ \t]+(\\d+)");

    if (_stream) {
        component.contigs.clear();
        component.gaps.clear();

        int state = eStart;
        std::string buf;
        boost::char_separator< char > sep(" ,\t");

        while (std::getline(_stream, buf)) {
            if (state == eStart) {
                boost::smatch what;
                if (boost::regex_match(buf, what, reg)) {
                    state = eId;
                    component.id = boost::lexical_cast< size_t >(what[1]);
                } else {
                    LOG4CXX_WARN(logger, boost::format("fa=>invalid line for not start with >: %s") % buf);
                    return false;
                }
            } else if (state == eId) {
                boost::tokenizer< boost::char_separator< char > > toker(buf, sep);
                BOOST_FOREACH(const std::string& id, toker)  {
                    component.contigs.push_back( boost::lexical_cast< size_t > (id));
                }
                state = eGap;
            } else if (state == eGap) {
                boost::tokenizer< boost::char_separator< char > > toker(buf, sep);
                BOOST_FOREACH(const std::string& gap, toker)  {
                    component.gaps.push_back( boost::lexical_cast< long > (gap));
                }

                component.gaps.push_back(0);//finall contige 
                if (component.contigs.size() == component.gaps.size()) {
                    state = eStart; 
                    return true;
                } else {
                    LOG4CXX_WARN(logger, boost::format("fa=>invalid line for id size != gap size: %s") % buf);
                    return false;
                }
            }
        }
    }
    return false;
}

std::ostream& operator<<(std::ostream& os, const Component& component) {
    os << "len:" << component._length << std::endl;
    for (size_t i=0; i<component.contigs.size(); ++i) {
        os << component.contigs[i] << "\t";
    }
    os << std::endl;
    for (size_t i=0; i<component.gaps.size(); ++i) {
        os << component.gaps[i] << "\t";
    }
    return os;
}

bool ReadComponents(std::istream& stream, const ContigList& contigs, size_t K, ComponentList& components) {
    LOG4CXX_DEBUG(logger, boost::format("read component from stream begin"));

    ComponentReader reader(stream);
    Component component;
    while (reader.read(component)) {
        component.length(K, contigs);
        components.push_back(component);
    }

    LOG4CXX_DEBUG(logger, boost::format("read component from stream end"));
    return true;
}

bool ReadComponents(const std::string& file, const ContigList& contigs, size_t K, ComponentList& components) {
    std::ifstream stream(file.c_str());
    return ReadComponents(stream, contigs, K, components);
}

