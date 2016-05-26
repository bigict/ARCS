#include "scaffolding.h"
#include "component.h"
#include "contigs.h"
#include "gapped_fragment_graph.h"
#include "kmer.h"
#include "kmer_tbl.h"
#include "kseq.h"
#include "pair_read.h"

#include <iostream>
#include <numeric>
#include <algorithm>
#include <tuple>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Scaffolding"));

template< size_t K >
class ConnectGraphBuilder {
public:
    ConnectGraphBuilder(size_t L, size_t insert_size, const PairReadList& pair_reads, const KmerMultiTable< K, KmerPosition >& hash_tbl, const ComponentList& components, bool do_reverse=true) : _K(L), _insert_size(insert_size), _pair_reads(pair_reads), _hash_tbl(hash_tbl), _components(components), _do_reverse(do_reverse), _read_len(-1) {
    }

    size_t build(size_t thread_num, GappedFragmentGraph* graph) {
        size_t num = 0;

        LOG4CXX_DEBUG(logger, boost::format("build graph begin"));
        if (thread_num <= 1) {
            num = callback(0, _pair_reads.size(), graph);
        } else {
            std::vector< size_t > pair_kmer_num(thread_num, 0);

            size_t block = _pair_reads.size() / thread_num;
            boost::thread_group group;
            for (size_t i = 0; i < thread_num; ++i) {
                size_t start = i * block, end = (i+1 == thread_num ? _pair_reads.size() : (i + 1) * block);
                group.create_thread(boost::bind((&ConnectGraphBuilder::callback), this, start, end, i, thread_num, graph, &pair_kmer_num[i]));
            }
            group.join_all();

            LOG4CXX_DEBUG(logger, boost::format("multi thread build graph end"));
            num = std::accumulate(pair_kmer_num.begin(), pair_kmer_num.end(), 0);
        }
        LOG4CXX_DEBUG(logger, boost::format("build graph end"));

        return num;
    }

    void callback(size_t start, size_t end, size_t i, size_t thread_num, GappedFragmentGraph * graph, size_t *pair_kmer_num) {
        LOG4CXX_DEBUG(logger, boost::format("start thread %d/%d start=%d end=%d") % i % thread_num % start % end);

        size_t num = callback(start, end, graph);
        if (pair_kmer_num != NULL) {
            *pair_kmer_num = num;
        }

        LOG4CXX_DEBUG(logger, boost::format("stop thread %d/%d start=%d end=%d") % i % thread_num % start % end);
    }

    size_t callback(size_t start, size_t end, GappedFragmentGraph * graph) {
        size_t pair_kmer_num = 0;

        for (size_t i = start; i < end; ++i) {
            const PairRead& pair_read = _pair_reads[i];
            if (pair_read.set1.size() < _K || pair_read.set2.size() < _K) {
                continue;
            }
            pair_kmer_num += addEdge(pair_read.set1, make_complement_dna(pair_read.set2), graph);
            if (_do_reverse) {
                pair_kmer_num += addEdge(pair_read.set2, make_complement_dna(pair_read.set1), graph);
            }
        }

        return pair_kmer_num;
    }

    size_t getReadLen() const {
        BOOST_ASSERT(_read_len != -1);
        return _read_len;
    }
private:

    size_t addEdge(const std::string& read1, const std::string& read2, GappedFragmentGraph* graph) {
        if(_read_len == -1) {
            _read_len = read1.length();
        }
        size_t num = 0;

        std::map< std::pair< size_t, size_t >, std::tuple<long, long, int> > find_component;
        int kmer_num = 0;
        bool first = false;
        std::set<size_t> find_contigs;
        for (size_t i = 0, j = _K; j <= read1.size() && j <= read2.size(); ++i,++j,++kmer_num) {
            Kmer< K > kmer1 = read1.substr(i, _K), kmer2 = read2.substr(i, _K);
            BOOST_ASSERT(_hash_tbl.count(kmer1) <= 1);
            BOOST_ASSERT(_hash_tbl.count(kmer2) <= 1);

            typename KmerMultiTable< K, KmerPosition >::const_iterator left = _hash_tbl.find(kmer1), right = _hash_tbl.find(kmer2);
            
            if(left != _hash_tbl.end() && right != _hash_tbl.end() && left->second.first == right->second.first) {
                if(find_contigs.find(left->second.first) == find_contigs.end()
                        && left->second.second <= 200 
                        && right->second.second <= 400) {
                    find_contigs.insert(left->second.first);
                    graph->node2PairReads[ left->second.first ] += 1;
                }
            }
            if (left != _hash_tbl.end() && right != _hash_tbl.end() && left->second.first != right->second.first) {
                ++num;

                LOG4CXX_TRACE(logger, boost::format("pair kmre1=%s %d %d") % kmer1 % left->second.first % left->second.second);
                LOG4CXX_TRACE(logger, boost::format("pair kmre2=%s %d %d") % kmer2 % right->second.first % right->second.second);
            
                //////////////////////////////////////////////////////
                //
                //              x
                //          C1  |                       y
                // |____________|____|             C2   |
                // |            |    |    |_____________|__|
                // |            |    | gap|             |
                // |            |    |----|             |
                // |            |        insert_size    |
                // |            |-----------------------|
                // |     distance         |
                // |----------------------|
                //
                // NOTE: overlap := -gap
                // 
                // ///////////////////////////////////////////////////
                long distance = _insert_size + left->second.second - right->second.second;
                size_t C1 = _components[left->second.first].length();
                long overlap = C1 - distance;
                if (overlap > 0 && overlap > _K) {
                    if (overlap > _insert_size)
                        continue;
                    distance = C1;
                }

                boost::unique_lock< boost::mutex > lock(_mtx);
                graph->addEdge(left->second.first, right->second.first, distance, 1, 0);
                std::pair<size_t, size_t> p = std::make_pair(left->second.first, right->second.first);
                if(find_component.find(p) == find_component.end()) {
                    find_component[p] = std::make_tuple(left->second.second, left->second.second, 1);
                } else {
                    std::get<0>(find_component[p]) = std::min(left->second.second, std::get<0>(find_component[p]));
                    std::get<1>(find_component[p]) = std::max(left->second.second, std::get<1>(find_component[p]));
                    std::get<2>(find_component[p]) += 1;
                }
            }
        }
        std::pair<size_t, size_t> pos(-1, -1);
        int mx = -1;
        for(auto it = find_component.begin(); it != find_component.end(); ++it) {
            if(std::get<2>(it->second) > mx && 
                    ( std::get<2>(it->second) > read1.length() - _K || 
                      std::get<2>(it->second) >= std::max(2.0, 0.8 * (std::get<1>(it->second) - std::get<0>(it->second))))
            ) {
                mx = std::get<2>(it->second);
                pos = it->first;
            }
        }
        if(mx != -1) {
            boost::unique_lock< boost::mutex > lock(_mtx);
            graph->addEdge(pos.first, pos.second, 0, 0, 1);
        }

        return num;
    }

    size_t _K;
    size_t _insert_size;
    const PairReadList& _pair_reads;
    const KmerMultiTable< K, KmerPosition >& _hash_tbl;
    const ComponentList& _components;
    bool _do_reverse;
    size_t _read_len;
    boost::mutex _mtx;
};

template< size_t K >
int _Scaffolding_run_(size_t L, const Properties& options, const Arguments& arguments) {
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
    if (!ReadComponents(options.get< std::string >("f"), contigs, components)) {
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

        InsertSizeEstimater< K > estimater(L, INSERT_SIZE, pair_reads, hash_tbl, options.find("S") == options.not_found());
        estimater.estimate(&INSERT_SIZE, &DELTA);
        if(INSERT_SIZE == 0) {
            if(options.find("L") == options.not_found()) {
                LOG4CXX_ERROR(logger, boost::format("specify INSERT_SIZE using -L or in the config file"));
                return 1;
            } else {
                INSERT_SIZE = options.get< size_t >("L");
            }
        }
        DELTA = std::max(DELTA, 1.0E-5);
    }
    LOG4CXX_INFO(logger, boost::format("INSERT_SIZE=[%d], DELTA=[%f],EDGE_CUTOFF=[%d]") % INSERT_SIZE % DELTA % EDGE_CUTOFF);
    //build graph
    size_t pair_kmer_cutoff = options.get< size_t >("r", 0); //no use "wangbing"
    size_t pair_read_cutoff = options.get< size_t >("R", 0);;
    double percent = options.get< double >("P", 0.0); // no use "wangbing"
    LOG4CXX_INFO(logger, boost::format("pair_kmer_cutoff=[%d], pair_read_cutoff=[%d],edge_link_percent=[%f]") % pair_kmer_cutoff % pair_read_cutoff % percent);

    GappedFragmentGraph g(L, pair_kmer_cutoff, pair_read_cutoff, percent, components.size(), GENOME_LEN);
    g.INSERT_SIZE = INSERT_SIZE;
    g.DELTA = DELTA;
    size_t read_len = -1;
    // build
    {
        KmerMultiTable< K, KmerPosition > hash_tbl;
        //BuildKmerTable_pairends< K >(L, options.get< size_t >("L", 180), EDGE_CUTOFF, contigs, components, hash_tbl);
        BuildKmerTable_pairends< K >(L, INSERT_SIZE, EDGE_CUTOFF, contigs, components, hash_tbl);

        ConnectGraphBuilder< K > builder(L, INSERT_SIZE, pair_reads, hash_tbl, components, options.find("S") == options.not_found());
        g.PAIR_KMER_NUM = builder.build(options.get< size_t >("p", 1), &g);
        read_len = builder.getReadLen();
    }
    //for_each(g.node2PairReads.begin(), g.node2PairReads.end(), [](const std::pair<size_t, int>& x){std::cout << x.second << std::endl;});
    LOG4CXX_INFO(logger, boost::format("pair_kmer_num=%d") % g.PAIR_KMER_NUM);

    // graph
    {
        boost::filesystem::ofstream stream(workdir / boost::filesystem::path(
                    boost::str(boost::format("contig_arc_graph_before_remove_ambigous_arcs_%d") % ITERATION)
            ));
        stream << g;
    }

    g.scoreAndRemoveNoise(components, read_len);

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
        g.outputLP(stream, components);
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
    //K = K - 1;

    // process
    if (0 < K && K <= 32) {
        r = _Scaffolding_run_< 32 >(K, options, arguments);
    } else if ( 32 < K && K <= 64) {
        r = _Scaffolding_run_< 64 >(K, options, arguments);
    } else if ( 64 < K && K <= 96) {
        r = _Scaffolding_run_< 96 >(K, options, arguments);
    } else if ( 96 < K && K <= 128) {
        r = _Scaffolding_run_< 128 >(K, options, arguments);
    } else if (128 < K && K <= 160) {
        r = _Scaffolding_run_< 160 >(K, options, arguments);
    } else if (160 < K && K <= 192) {
        r = _Scaffolding_run_< 192 >(K, options, arguments);
    }

    LOG4CXX_DEBUG(logger, "scaffolding end");
    return r;
}

Scaffolding Scaffolding::_runner;

Scaffolding::Scaffolding() : Runner("d:K:C:f:e:1:2:L:P:p:i:r:R:c:Sh") {
    RUNNER_INSTALL("scaffold", this, "generate ordered sets of contigs using distance estimates");
}

int Scaffolding::printHelps() const {
    std::cout << "usage: arcs scaffold [arguments]" << std::endl;
    std::cout << std::endl;
    std::cout << "\t-d[=<workdir>] word dir, default ." << std::endl;
    std::cout << "\t-K[=<number>]  kmer size, default 31" << std::endl;
    std::cout << "\t-C[=<file>]    contig file <cdbg_copy_number.fa>" << std::endl;
    std::cout << "\t-f[=<file>]    component file <component_\%d>" << std::endl;
    std::cout << "\t-e[=<number>]  edge length cutoff, ignore component whose length < edge length cutoff, default number of -L" << std::endl;
    std::cout << "\t-1[=<file>]    pair reads file 1" << std::endl;
    std::cout << "\t-2[=<file>]    pair reads file 2" << std::endl;
    std::cout << "\t-L[=<number>]  insert size, default 180 <arcs estimates a new insert size based on -L>" << std::endl;
    //std::cout << "\t-P" << std::endl;
    std::cout << "\t-p[=<number>]  thread number, default 1" << std::endl;
    std::cout << "\t-i[=<number>]  number of iterate time of scaffolding" << std::endl;
    //std::cout << "\t-r" << std::endl;
    std::cout << "\t-R[=<number>]  pair read number threshold to build edges between components, default 0" << std::endl;
    std::cout << "\t-c[=<file>]    log config file, default ./log4cxx.properties" << std::endl;
    std::cout << "\t-S             do not reverse original reads" << std::endl;
    std::cout << "\t-h             help" << std::endl;
    std::cout << std::endl;
    return 256;
}

int Scaffolding::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
