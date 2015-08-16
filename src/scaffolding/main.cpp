#include "graph.h"
#include "kmer_tbl.h"
#include "pair_read.h"

#include <numeric>
#include <fstream>
#include <iostream>
#include <vector>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

typedef boost::property_tree::ptree Properties;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Contiging.main"));


class Scaffoding {
public:
    int run(const Properties& options) {
        if (checkOptions(options) != 0) {
            return 1;
        }

        // read contig file
        LOG4CXX_INFO(logger, boost::format("read contig file begin"));
        std::string contig_file_name = options.get<std::string>("C");
        size_t K = options.get< size_t >("K");
        LOG4CXX_INFO(logger, boost::format("Kmer size:K =  %d") % K);
        std::ifstream contig_in(contig_file_name.c_str());
        ContigSet contigs(contig_in, K);
        LOG4CXX_INFO(logger, boost::format("all contig length GENOME_LEN: %d") % contigs.GENOME_LEN);
        LOG4CXX_INFO(logger, boost::format("read contig file end"));

        // read component file
        LOG4CXX_INFO(logger, boost::format("read component file begin"));
        std::vector< Component > components;
        std::string component_file_name = options.get<std::string>("f");
        std::ifstream component_in(component_file_name.c_str());
        ComponentReader reader(component_in);
        Component component(K);
        while(reader.read(component)){
            component.initializeLen(contigs);
            components.push_back(component);
        }
        LOG4CXX_INFO(logger, boost::format("Components size : %d") % components.size());
        LOG4CXX_INFO(logger, boost::format("read component file end"));

        // build kmer table
        LOG4CXX_INFO(logger, boost::format("build kmer table for insert size begin"));
        size_t INSERT_SIZE = 180;
        if (options.find("L") != options.not_found()) {
            INSERT_SIZE = options.get<size_t>("L");
        }
        size_t EDGE_CUTOFF = K;
        if (options.find("e") != options.not_found()) {
            EDGE_CUTOFF = options.get<size_t>("e");
        }
        KmerTable table_for_insert_size(K);
        size_t count = 0;
        BOOST_FOREACH(Component& component, components) {
            LOG4CXX_TRACE(logger, boost::format("Component No = %d, size = %d, length = %d") % count % component._contig_id.size() %  contigs._contigs[component._contig_id[0]].contig.length());
            if (component._contig_id.size() == 0 || (component.getLen() < 2 * INSERT_SIZE)) {
                count ++;
                continue;
            }
            component.produceKmerForInsertSize(contigs, table_for_insert_size, count);
            count ++;
        }
        LOG4CXX_INFO(logger, boost::format("FOR INSERT SIZE kmer number = %d, EDGE_CUTOFF = %d") % table_for_insert_size.size() % EDGE_CUTOFF);
        LOG4CXX_INFO(logger, boost::format("build kmer table for insert size end"));

        LOG4CXX_INFO(logger, boost::format("build kmer table for pair read begin"));
        KmerTable table_for_pair_read(K);
        count = 0;
        BOOST_FOREACH(Component& component, components) {
            if (component._contig_id.size() == 0 || (component._contig_id.size() == 1  && contigs._contigs[component._contig_id[0]].contig.length() < EDGE_CUTOFF)) {
            //if (component._contig_id.size() == 0 || (component.getLen() < EDGE_CUTOFF)) {
                count ++;
                continue;
            }
            component.produceKmerForPairRead(contigs, table_for_pair_read, count, INSERT_SIZE);
            count ++;
        }
        LOG4CXX_INFO(logger, boost::format("FOR PAIR READ kmer number = %d, EDGE_CUTOFF = %d") % table_for_pair_read.size() % EDGE_CUTOFF);
        LOG4CXX_INFO(logger, boost::format("build kmer table for pair read end"));


        LOG4CXX_INFO(logger, boost::format("read pair read begin"));
        std::string read_file1 = options.get<std::string>("1");
        std::string read_file2 = options.get<std::string>("2");
        std::ifstream r1(read_file1.c_str());
        std::ifstream r2(read_file2.c_str());

        PairReadSet pair_read(r1, r2, K, INSERT_SIZE);
        LOG4CXX_INFO(logger, boost::format("pair read num = %d") % pair_read.size());
        LOG4CXX_INFO(logger, boost::format("read pair read end"));

        //build graph
        LOG4CXX_INFO(logger, boost::format("build graph begin"));
        size_t pair_kmer_cutoff = 0;
        size_t pair_read_cutoff = 0;
        double percent = 0.0;
        if (options.find("r") != options.not_found()) {
            pair_kmer_cutoff = options.get<size_t>("r");
        }
        if (options.find("R") != options.not_found()) {
            pair_read_cutoff = options.get<size_t>("R");
        }
        if (options.find("P") != options.not_found()) {
            percent = options.get<double>("P");
        }
        pair_read.estimateInsertSize(table_for_insert_size);
        Graph g(K, pair_kmer_cutoff, pair_read_cutoff, percent, components.size(), contigs.GENOME_LEN);
        pair_read.buildConnectGraph(g, table_for_pair_read, components);
        g.setPairKmerNumAndInsertSizeAndDelta(pair_read.PAIR_KMER_NUM, pair_read.INSERT_SIZE, pair_read.DELTA);
        LOG4CXX_INFO(logger, boost::format("pair_kmer_num=%d") % pair_read.PAIR_KMER_NUM);
        g.scoreAndRemoveNoise(components);
        //g.outputLP(std::cout);
        std::cout << g << std::endl;
        LOG4CXX_INFO(logger, boost::format("build graph end"));
    }
    private:
        int checkOptions(const Properties& options);
        int printHelps() const;
};

int Scaffoding::checkOptions(const Properties& options) {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }

    std::vector< std::string > necessary = boost::assign::list_of("d")("K")("C")("f")("1")("2");
    BOOST_FOREACH(const std::string& c, necessary) {
        if (options.find(c) == options.not_found()) {
            std::cerr << boost::format("The argument -%s is necessary! type -h for more help") % c << std::endl;
            return 1;
        }
    }

    std::vector< std::string > intopt = boost::assign::list_of("K");
    BOOST_FOREACH(const std::string& c, intopt) {
        try {
            size_t K = options.get< size_t >(c);
        } catch (std::exception& e) {
            std::cerr << boost::format("The argument -%s should be a integer. type -h for more help") % c << std::endl;
            return 1;
        }
    }

    return 0;
}

int Scaffoding::printHelps() const {
    std::cerr << "Scaffording - version v0.9" << std::endl
        << "Usage:" << std::endl
        << "Arguments:" << std::endl
        << "-d:     select your work directory." << std::endl
        << "-K:     select your kmer size." << std::endl
        << "-C:     select a contig file" << std::endl
        << "-f:     select a component file" << std::endl
        << "-e:     contig length cutoff" << std::endl
        << "-1:     read file 1" << std::endl
        << "-2:     read file 2" << std::endl
        << "-L:     insert size" << std::endl
        << "-P:     link qulity percent" << std::endl
        << "-p:     select the number of the threads you want to use. The argument--p must be a figure more than zero." << std::endl
        << "-i:     iterator time" << std::endl
        << "-r:     pair kmer cutoff" << std::endl
        << "-R:     pair read cutoff" << std::endl
        << "-c:     select a log config file." << std::endl
        << "-h:     help" << std::endl;

    return 1;
}

int main(int argc, char* argv[]) {
    Properties options;
    {
        // command line options
        Properties cmd;
        std::string opt_string("d:K:C:f:e:1:2:L:P:p:i:r:R:c:h");
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

        // merge options
        for (Properties::const_iterator it = cmd.begin(); it != cmd.end(); it++){
            options.put(it->first,it->second.data());
        }
    }

    // build contigs.
    Scaffoding c;
    return c.run(options);
}

