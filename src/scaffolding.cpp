#include "scaffolding.h"
#include "component.h"
#include "contigs.h"
#include "gapped_fragment_graph.h"
#include "kmer.h"
#include "kmer_tbl.h"
#include "kseq.h"
#include "pair_read.h"

#include <iostream>

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Scaffolding"));
Scaffolding Scaffolding::_runner;

template< size_t K >
class ConnectGraphBuilder {
public:
    ConnectGraphBuilder(size_t L, size_t insert_size, const PairReadList& pair_reads, const KmerTable< K, KmerPosition >& hash_tbl, const ComponentList& components) : _K(L), _insert_size(insert_size), _pair_reads(pair_reads), _hash_tbl(hash_tbl), _components(components) {
    }

    size_t build(GappedFragmentGraph* graph) const {
        size_t pair_kmer_num = 0;

        LOG4CXX_DEBUG(logger, boost::format("build graph begin"));
        BOOST_FOREACH(const PairRead& pair_read, _pair_reads) {
            if (pair_read.set1.size() < _K || pair_read.set2.size() < _K) {
                continue;
            }
            pair_kmer_num += addEdge(pair_read.set1, make_complement_dna(pair_read.set2), graph);
            pair_kmer_num += addEdge(pair_read.set2, make_complement_dna(pair_read.set1), graph);
        }
        LOG4CXX_DEBUG(logger, boost::format("build graph end"));

        return pair_kmer_num;
    }
private:
    struct PairKmerCmp {
        bool operator()(const std::pair<size_t, size_t>& l, const std::pair<size_t, size_t>& r) {
            if (l.first != r.first) {
                return l.first < r.first;
            }
            return l.second < r.second;
        }
    };

    size_t addEdge(const std::string& read1, const std::string& read2, GappedFragmentGraph* graph) const {
        size_t num = 0;

        std::set< std::pair< size_t, size_t >, PairKmerCmp > find_component;
        for (size_t i = 0, j = _K; j <= read1.size() && j <= read2.size(); ++i,++j) {
            Kmer< K > kmer1 = read1.substr(i, _K), kmer2 = read2.substr(i, _K);

            typename KmerTable< K, KmerPosition >::const_iterator left = _hash_tbl.find(kmer1), right = _hash_tbl.find(kmer2);
            
            if (left != _hash_tbl.end() && right != _hash_tbl.end() && left->second.first != right->second.first) {
                ++num;

                LOG4CXX_TRACE(logger, boost::format("pair kmre1=%s %d %d") % kmer1 % left->second.first % left->second.second);
                LOG4CXX_TRACE(logger, boost::format("pair kmre2=%s %d %d") % kmer2 % right->second.first % right->second.second);
            
                long distance = _insert_size + left->second.second - right->second.second;
                long left_len = _components[left->second.first].length();
                long overlap = left_len - distance - _K + 1;
                if (overlap > 0) {
                    if (overlap > _insert_size)
                        continue;
                    distance = left_len - _K + 1; //can change
                }
                if (find_component.find(std::make_pair(left->second.first, right->second.first)) != find_component.end()) {
                    graph->addEdge(left->second.first, right->second.first, distance, 1, 0);
                } else {
                    find_component.insert( std::make_pair(left->second.first, right->second.first));
                    graph->addEdge(left->second.first, right->second.first, distance, 1, 1);
                }
            }
        }

        return num;
    }

    size_t _K;
    size_t _insert_size;
    const PairReadList& _pair_reads;
    const KmerTable< K, KmerPosition >& _hash_tbl;
    const ComponentList& _components;
};

template< size_t K >
int _Scaffolding_run_(size_t L, const Properties& options) {
    // check root dir
    boost::filesystem::path workdir(options.get< std::string >("d", "."));
    if (!boost::filesystem::exists(workdir) && !boost::filesystem::create_directory(workdir)) {
        LOG4CXX_ERROR(logger, boost::format("failed to create directory: %s") % workdir);
        return 1;
    }

    Kmer< K >::length(L); // IMPORTANT: set kmer length

    size_t ITERATION = options.get< size_t >("i");

    // read contig file
    ContigList contigs;
    if (!ReadContigs(options.get< std::string >("C"), contigs)) {
        LOG4CXX_ERROR(logger, boost::format("faild to read contigs from file: %s") % options.get< std::string >("C"));
        return 1;
    }
    size_t GENOME_LEN = GenomeLength(contigs, L);
    LOG4CXX_INFO(logger, boost::format("GENOME_LEN=%d") % GENOME_LEN);

    // read component file
    ComponentList components;
    if (!ReadComponents(options.get< std::string >("f"), contigs, L, components)) {
        LOG4CXX_ERROR(logger, boost::format("faild to read components from file: %s") % options.get< std::string >("f"));
        return 1;
    }
    LOG4CXX_INFO(logger, boost::format("Components size=%d") % components.size());

    size_t INSERT_SIZE = options.get< size_t >("L", 180);;
    double DELTA = 0;
    size_t EDGE_CUTOFF = options.get< size_t >("e", L);
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
        KmerTable< K, KmerPosition > hash_tbl;
        BuildKmerTable_insertsize< K >(L, INSERT_SIZE, contigs, components, hash_tbl);

        InsertSizeEstimater< K > estimater(K, INSERT_SIZE, pair_reads, hash_tbl);
        estimater.estimate(&INSERT_SIZE, &DELTA);
    }
    LOG4CXX_INFO(logger, boost::format("INSERT_SIZE=[%d], DELTA=[%f]") % INSERT_SIZE % DELTA);
    //build graph
    size_t pair_kmer_cutoff = options.get< size_t >("r", 0);
    size_t pair_read_cutoff = options.get< size_t >("R", 0);;
    double percent = options.get< double >("P", 0.0);

    GappedFragmentGraph g(L, pair_kmer_cutoff, pair_read_cutoff, percent, components.size(), GENOME_LEN);
    // build
    {
        KmerTable< K, KmerPosition > hash_tbl;
        BuildKmerTable_pairends< K >(L, options.get< size_t >("L", 180), EDGE_CUTOFF, contigs, components, hash_tbl);

        ConnectGraphBuilder< K > builder(L, INSERT_SIZE, pair_reads, hash_tbl, components);
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

int Scaffolding::run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    LOG4CXX_DEBUG(logger, "scaffolding begin");

    size_t K = options.get< size_t >("K", 31);

    // process
    if (0 < K && K <= 32) {
        r = _Scaffolding_run_< 32 >(K, options);
    } else if (32 < K && K <= 64) {
        r = _Scaffolding_run_< 64 >(K, options);
    } else if (64 < K && K <= 96) {
        r = _Scaffolding_run_< 64 >(K, options);
    } else if (96 < K && K <= 128) {
        r = _Scaffolding_run_< 128 >(K, options);
    } else if (96 < K && K <= 160) {
        r = _Scaffolding_run_< 160 >(K, options);
    } else if (160 < K && K <= 192) {
        r = _Scaffolding_run_< 192 >(K, options);
    }

    LOG4CXX_DEBUG(logger, "scaffolding end");
    return r;
}

Scaffolding::Scaffolding() : Runner("d:K:C:f:e:1:2:L:P:p:i:r:R:c:h") {
    RUNNER_INSTALL("scaffolding", this, "scaffold");
}

int Scaffolding::printHelps() const {
    std::cout << "arcs scaffolding -K [kmer] -i [input] -d [workdir]" << std::endl;
    return 256;
}

int Scaffolding::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
