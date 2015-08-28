#include "contigs.h"

#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>


#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scanfording.contig"));

ContigSet::ContigSet(std::istream& stream, size_t K) : _stream(stream), _k(K) {
    if (!stream) {
        LOG4CXX_INFO(logger,boost::format("can not find contig file"));
    }
    GENOME_LEN = 0;
    Contig contig;
    while (read(contig)) {
        GENOME_LEN += contig.contig.size() - _k + 1;
        _contigs.push_back(contig);  
    }
}

bool ContigSet::read(Contig& contig) {
    enum {
        eCopy,
        eContig, 
    };

    if (_stream) {
        int state = eCopy;
        std::string buf;
        boost::char_separator<char> fn(" ,\t");

        while (std::getline(_stream, buf)) {
            if (state == eCopy) {
                if (boost::algorithm::starts_with(buf, ">")) {
                    boost::tokenizer< boost::char_separator<char> > toker(buf, fn);
                //    if(toker.size()==2)
                //    {
                    boost::tokenizer< boost::char_separator<char> >::iterator  iter = toker.begin();  
                    contig.copy_num = boost::lexical_cast<int> (*(++iter) );
                    state = eContig;  
                //    }else{
                //        LOG4CXX_WARN(logger, boost::format("fa=>invalid line for copy_num: %s") % buf);
                //        return false;
                //    }

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

std::ostream& operator<<(std::ostream& os, const ContigSet& contigs) {
    size_t index = 0;
    BOOST_FOREACH(const Contig& contig, contigs._contigs) {
        os << index << "\t" << contig.contig << "\t" << contig.copy_num << "\t" << contig.contig.length() << std::endl;
        ++ index; 
    }
    return os;
} 

