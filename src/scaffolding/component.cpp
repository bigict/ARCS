#include "component.h"
#include "kmer_tbl.h"

#include <fstream>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffolding.component"));

size_t Component::produceKmerForInsertSize(const ContigList& contigs, KmerTable& tbl, size_t component_no) {
    long index = 0;
    size_t tbl_size = 0; 
    for (size_t j=0; j<_contig_id.size(); j++) {
        size_t con_id = _contig_id[j];
        for (size_t i = 0; i < contigs[con_id].seq.length() - _K + 1; i++) {
            Kmer kmer(contigs[con_id].seq, i, i + _K);
            tbl.addKmer(kmer, std::make_pair(component_no, index)); 
            index ++;
            tbl_size ++;
        }
        index += _gap[j]; 
    }
    return tbl_size;
}

size_t Component::produceKmerForPairRead(const ContigList& contigs, KmerTable& tbl, size_t component_no, size_t INSERT_SIZE) {
    long index = 0; 
    size_t tbl_size = 0;
    int cutoff = 2*INSERT_SIZE;
    if (_contig_id.size() == 1) {
        size_t con_id = _contig_id[0];
        for (size_t i = 0; i < contigs[con_id].seq.length() - _K + 1; i++) {
            if (index <= cutoff || index >= _length - 1 - _K + 1 - cutoff) {
                Kmer kmer(contigs[con_id].seq, i, i + _K);
                tbl.addKmer(kmer, std::make_pair(component_no, index)); 
                tbl_size ++; 
            }
            index ++;
        }
    } else {
        for (size_t j=0; j<_contig_id.size(); j++) {
            size_t con_id = _contig_id[j];
            for (size_t i = 0; i < contigs[con_id].seq.length() - _K + 1; i++) {
                if (index <= cutoff || index >= _length - 1 - cutoff) {
                    Kmer kmer(contigs[con_id].seq, i, i + _K);
                    tbl.addKmer(kmer, std::make_pair(component_no, index)); 
                    tbl_size ++;
                    std::cerr << component_no <<" "<< index << std::endl;
                }
                index ++;
            }
            index += _gap[j]; 
        }
    }
    return tbl_size;
}

void Component::length(const ContigList& contigs) {
    _length = 0;
    if (!contigs.empty()) {
        _length = contigs[_contig_id[0]].seq.length() + 1;//can change
        for (size_t i = 1; i < _contig_id.size(); ++i) {
            _length += _gap[i - 1];
            BOOST_ASSERT(contigs[_contig_id[i]].seq.size() >= _K - 1);
            _length += contigs[_contig_id[i]].seq.size() - _K + 1;
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
        component._contig_id.clear();
        component._gap.clear();

        int state = eStart;
        std::string buf;
        boost::char_separator< char > sep(" ,\t");

        while (std::getline(_stream, buf)) {
            if (state == eStart) {
                boost::smatch what;
                if (boost::regex_match(buf, what, reg)) {
                    state = eId;
                    // component.id = boost::lexical_cast< size_t >(what[1]);
                } else {
                    LOG4CXX_WARN(logger, boost::format("fa=>invalid line for not start with >: %s") % buf);
                    return false;
                }
            } else if (state == eId) {
                boost::tokenizer< boost::char_separator< char > > toker(buf, sep);
                BOOST_FOREACH(const std::string& id, toker)  {
                    component._contig_id.push_back( boost::lexical_cast< size_t > (id));
                }
                state = eGap;
            } else if (state == eGap) {
                boost::tokenizer< boost::char_separator< char > > toker(buf, sep);
                BOOST_FOREACH(const std::string& gap, toker)  {
                    component._gap.push_back( boost::lexical_cast< long > (gap));
                }

                component._gap.push_back(0);//finall contige 
                if (component._contig_id.size() == component._gap.size()) {
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
    for (size_t i=0; i<component._contig_id.size(); ++i) {
        os << component._contig_id[i] << "\t";
    }
    os << std::endl;
    for (size_t i=0; i<component._gap.size(); ++i) {
        os << component._gap[i] << "\t";
    }
    return os;
}

bool ReadComponents(std::istream& stream, const ContigList& contigs, size_t K, ComponentList& components) {
    LOG4CXX_DEBUG(logger, boost::format("read component from stream begin"));

    ComponentReader reader(stream);
    Component component(K);
    while (reader.read(component)) {
        component.length(contigs);
        components.push_back(component);
    }

    LOG4CXX_DEBUG(logger, boost::format("read component from stream end"));
    return true;
}

bool ReadComponents(const std::string& file, const ContigList& contigs, size_t K, ComponentList& components) {
    std::ifstream stream(file.c_str());
    return ReadComponents(stream, contigs, K, components);
}

