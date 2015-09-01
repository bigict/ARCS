#include "graph.h"
#include "kmer_tbl.h"
#include "pair_read.h"

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

typedef boost::property_tree::ptree Properties;

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffolding.main"));


class Scaffoding {
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
        size_t ITERATION = options.get< size_t >("i");

        LOG4CXX_INFO(logger, boost::format("K=%d") % K);

        // read contig file
        ContigList contigs;
        if (!ReadContigs(options.get< std::string >("C"), contigs)) {
            LOG4CXX_ERROR(logger, boost::format("faild to read contigs from file: %s") % options.get< std::string >("C"));
            return 1;
        }
        size_t GENOME_LEN = GenomeLength(contigs, K);
        LOG4CXX_INFO(logger, boost::format("GENOME_LEN=%d") % GENOME_LEN);

        // read component file
        ComponentList components;
        if (!ReadComponents(options.get< std::string >("f"), contigs, K, components)) {
            LOG4CXX_ERROR(logger, boost::format("faild to read components from file: %s") % options.get< std::string >("f"));
            return 1;
        }
        LOG4CXX_INFO(logger, boost::format("Components size=%d") % components.size());

        size_t INSERT_SIZE = options.get< size_t >("L", 180);;
        double DELTA = 0;
        size_t EDGE_CUTOFF = options.get< size_t >("e", K);
        size_t PAIR_KMER_NUM = 0;

        // read pair reads file
        PairReadList pair_reads;
        if (!ReadPairReads(options.get< std::string >("1"), options.get< std::string >("2"), pair_reads)) {
            LOG4CXX_ERROR(logger, boost::format("faild to read components from file: <%s,%s>") % options.get< std::string >("1") % options.get< std::string >("2"));
            return 1;
        }
        LOG4CXX_INFO(logger, boost::format("pair read num = %d") % pair_reads.size());

        // estimate insert size
        {
            KmerList hash_tbl;
            BuildKmerTable_insertsize(K, INSERT_SIZE, contigs, components, hash_tbl);

            InsertSizeEstimater estimater(K, INSERT_SIZE, pair_reads, hash_tbl);
            estimater.estimate(&INSERT_SIZE, &DELTA);
        }
        LOG4CXX_INFO(logger, boost::format("INSERT_SIZE=[%d], DELTA=[%f]") % INSERT_SIZE % DELTA);

        //build graph
        size_t pair_kmer_cutoff = options.get< size_t >("r", 0);
        size_t pair_read_cutoff = options.get< size_t >("R", 0);;
        double percent = options.get< double >("P", 0.0);

        Graph g(K, pair_kmer_cutoff, pair_read_cutoff, percent, components.size(), GENOME_LEN);

        // build
        {
            KmerList hash_tbl;
            BuildKmerTable_pairends(K, options.get< size_t >("L", 180), EDGE_CUTOFF, contigs, components, hash_tbl);

            ConnectGraphBuilder builder(K, INSERT_SIZE, pair_reads, hash_tbl, components);
            PAIR_KMER_NUM = builder.build(&g);
        }
        LOG4CXX_INFO(logger, boost::format("pair_kmer_num=%d") % PAIR_KMER_NUM);

        g.setPairKmerNumAndInsertSizeAndDelta(PAIR_KMER_NUM, INSERT_SIZE, DELTA);

        // graph
        {
            boost::filesystem::ofstream stream(workdir / boost::filesystem::path(
                        boost::str(boost::format("contig_arc_graph_before_remove_ambigous_arcs_%d") % ITERATION)
                ));
            stream << g;
        }

        g.scoreAndRemoveNoise(components);

        // graph
        {
            boost::filesystem::ofstream stream(workdir / boost::filesystem::path(
                        boost::str(boost::format("contig_arc_graph_after_remove_ambigous_arcs_%d") % ITERATION)
                ));
            stream << g;
        }

        // parameters
        {
            boost::filesystem::ofstream stream(workdir / boost::filesystem::path(
                        boost::str(boost::format("scaffold_parameter_%d") % ITERATION)
                ));
            stream << boost::format("K=%d") % K << std::endl;
            stream << boost::format("EDGE_CLUSTER_NUM=%d") % components.size() << std::endl;
            stream << boost::format("GENOME_LEN=%d") % g.GENOME_LEN << std::endl;
            stream << boost::format("INSERT_SIZE=%d") % g.INSERT_SIZE << std::endl;
            stream << boost::format("DELTA=%d") % g.DELTA << std::endl;
        }
        // LP
        {
            boost::filesystem::ofstream stream(workdir / boost::filesystem::path(
                        boost::str(boost::format("position_lp_%d.math") % ITERATION)
                ));
            g.outputLP(stream);
        }
        // components
        {
            boost::filesystem::ofstream stream(workdir / boost::filesystem::path(
                        boost::str(boost::format("edge_cluster_len_%d") % ITERATION)
                ));
            BOOST_FOREACH(Component& component, components) {
                stream << component.length() << std::endl;
            }
        }

        return 0;
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

    std::vector< std::string > intopt = boost::assign::list_of("KrReL");
    BOOST_FOREACH(const std::string& c, intopt) {
        if (options.find(c) == options.not_found()) continue;
        try {
            options.get< size_t >(c);
        } catch (std::exception& e) {
            std::cerr << boost::format("The argument -%s should be a integer. type -h for more help") % c << std::endl;
            return 1;
        }
    }

    std::vector< std::string > doubleopt = boost::assign::list_of("P");
    BOOST_FOREACH(const std::string& c, doubleopt) {
        if (options.find(c) == options.not_found()) continue;
        try {
            options.get< double >(c);
        } catch (std::exception& e) {
            std::cerr << boost::format("The argument -%s should be a double. type -h for more help") % c << std::endl;
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

