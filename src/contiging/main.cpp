#include <numeric>
#include <fstream>
#include <iostream>
#include <vector>

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

#include "debruijn_graph.h"
#include "kmer.h"
#include "kmer_tbl.h"
#include "kseq.h"

typedef boost::property_tree::ptree Properties;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("contiging.main"));

class ReadQuality {
public:
    ReadQuality(size_t max_count = -1) : _max_count(max_count) {
    }

    size_t quality(const std::vector< std::string >& filelist, size_t* count) const {
        size_t q = 0, n = 0;

        BOOST_FOREACH(const std::string& file, filelist) {
            std::ifstream stream (file.c_str());
            q += quality(stream, &n);
        }

        if (count != NULL) {
            *count += n;
        }
        return q;
    }
    size_t quality(std::istream& stream, size_t* count) const {
        size_t q = 0, n = 0;

        if (stream) {
            DNASeqReader reader(stream);
            DNASeq read;
            
            while (reader.read(read) && n < _max_count) {
                q += std::accumulate(read.quality.begin(), read.quality.end(), 0);
                n += read.seq.length();
            }
        }

        if (count != NULL) {
            *count += n;
        }
        return q;
    }

    double threshold(const std::vector< std::string >& filelist, double* min_threshold) const {
        size_t n = 0;
        size_t q = quality(filelist, &n);
        return threshold(q, n, min_threshold);
    }

    double threshold(std::istream& stream, double* min_threshold) const {
        size_t n = 0;
        size_t q = quality(stream, &n);
        return threshold(q, n, min_threshold);
    }
    
    double threshold(size_t quality, size_t count, double* min_threshold) const {
        if (count > 0) {
            if (min_threshold != NULL) {
                *min_threshold = (quality / count) * 0.8;
            }
            return (quality / count) * 0.9;
        }
        return 0;
    }

private:
    size_t _max_count;
};

class Contiging {
public:
    int run(const Properties& options) {
        if (checkOptions(options) != 0) {
            return 1;
        }

        boost::filesystem::path workdir(options.get< std::string >("d"));
        if (!boost::filesystem::exists(workdir) && !boost::filesystem::create_directory(workdir)) {
            LOG4CXX_ERROR(logger, boost::format("failed to create directory: %s") % workdir);
            return 1;
        }

        size_t K = options.get< size_t >("K");
        std::vector< std::string > filelist = boost::assign::list_of(options.get< std::string >("q1"))(options.get< std::string >("q2"));

        double avg_quality = 0, min_quality = 0, percent = 1.0;
        if (options.find("E") != options.not_found()) {
            ReadQuality statistics;
            avg_quality = statistics.threshold(filelist, &min_quality);
            percent = 0.95;
        }
        LOG4CXX_INFO(logger, boost::format("avg_quality=%lf min_quality=%lf") % avg_quality % min_quality);

        KmerTable tbl(K, avg_quality, min_quality, percent, options.get< size_t >("READ_LENGTH_CUTOFF", 100), options.find("S") == options.not_found());

        if (!tbl.read(filelist)) {
            LOG4CXX_ERROR(logger, boost::format("faild to open file: %s") % options.get< std::string >("q1"));
            return -1;
        }

        DeBruijnGraph graph(tbl);
        graph.compact();
        // edge
        {
            boost::filesystem::ofstream stream(workdir / boost::filesystem::path("condensed_de_bruijn_graph_before_trimming.data"));
            stream << graph;
        }

        graph.removeNoise();
        graph.removeNoise();

        graph.compact();
        // edge
        {
            boost::filesystem::ofstream stream(workdir / boost::filesystem::path("condensed_de_bruijn_graph_after_trimming.data"));
            stream << graph;
        }

        return 0;
    }
private:
    int checkOptions(const Properties& options);
    int printHelps() const;
};

int Contiging::checkOptions(const Properties& options) {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }

    std::vector< std::string > necessary = boost::assign::list_of("d")("s")("p")("K");
    BOOST_FOREACH(const std::string& c, necessary) {
        if (options.find(c) == options.not_found()) {
            std::cerr << boost::format("The argument -%s is necessary! type -h for more help") % c << std::endl;
            return 1;
        }
    }
    
    std::vector< std::string > intopt = boost::assign::list_of("K")("READ_LENGTH_CUTOFF");
    BOOST_FOREACH(const std::string& c, intopt) {
        try {
            if (options.find(c) != options.not_found()) {
                size_t K = options.get< size_t >(c);
            }
        } catch (std::exception& e) {
            std::cerr << boost::format("The argument -%s should be a integer. type -h for more help") % c << std::endl;
            return 1;
        }
    }
    
    return 0;
}

int Contiging::printHelps() const {
     std::cerr << "ARCS - version v0.9" << std::endl
         << "Usage:" << std::endl
         << "Arguments:" << std::endl
         << "-d:     select your work directory." << std::endl
         << "-s:     select your configure file name." << std::endl
         << "-E:     select whether you want to filter out bad reads." << std::endl
         << "-p:     select the number of the threads you want to use. The argument--p must be a figure more than zero." << std::endl
         << "-K:     select your kmer size." << std::endl
         << "-S:     select whether you want to read pair ends." << std::endl
         << "-c:     select a log config file." << std::endl
         << "-h:     help" << std::endl;

    return 1;
}

int main(int argc, char* argv[]) {
    Properties options;
    {
        // command line options
        Properties cmd;
        std::string opt_string("s:d:K:ESp:c:h");
        int opt = -1;
        while ((opt = getopt(argc, argv, opt_string.c_str())) != -1) {
            std::string key(1, (char)opt);
            if (optarg == NULL) {
                cmd.put(key, NULL);
            } else {
                cmd.put(key, optarg);
            }
        }

        // config log4cxx.
        if (cmd.find("c") != cmd.not_found()) {
            log4cxx::PropertyConfigurator::configure(cmd.get< std::string >("c"));
        } else {
            log4cxx::BasicConfigurator::configure();
        }
        
        // load ini options
        if (cmd.find("s") != cmd.not_found()) {
            try {
                boost::property_tree::read_ini(cmd.get< std::string >("s"), options);
            } catch (const boost::property_tree::ini_parser_error& e) {
                LOG4CXX_ERROR(logger, boost::format("load test.cfg failed(%s).") % e.what());
                return 1;
            }
        }

        // merge options
        for (Properties::const_iterator it = cmd.begin(); it != cmd.end(); it++){
            options.put(it->first,it->second.data());
        }
    }

    // build contigs.
    Contiging c;
    return c.run(options);
}
