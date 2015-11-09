bin_PROGRAMS=arcs
arcs_SOURCES=arff.cpp \
            arff.h \
            component.cpp \
            component.h \
            condensed_debruijn_graph.cpp \
            condensed_debruijn_graph.h \
            condensed_debruijn_graph_reader.cpp \
            condensed_debruijn_graph_reader.h \
            contiging.cpp \
            contiging.h \
            contigs.cpp \
            contigs.h \
            constant.cpp \
            constant.h \
            copy_num_estimator.cpp \
            copy_num_estimator.h \
            debruijn_graph.h \
            gap_filler.cpp \
            gap_filler.h \
            gap_filling.cpp \
            gap_filling.h \
            gapped_fragment_graph.cpp \
            gapped_fragment_graph.h \
            contigs.h \
            kmer.h \
            kmer_dataset.h \
            kmer_tbl.h \
            kseq.cpp \
            kseq.h \
            lp_solver.cpp \
            lp_solver.h \
            pair_read.cpp \
            pair_read.h \
            preprocess.cpp \
            preprocess.h \
            runner.h \
            remove_repeats.cpp \
            remove_repeats.h \
            scaffolding.cpp \
            scaffolding.h \
            uniq_edge_graph.cpp \
            uniq_edge_graph.h \
            utils.h \
            main.cpp
arcs_CXXFLAGS=${LOG4CXX_CFLAGS}
arcs_LDFLAGS=${LOG4CXX_LIBS} \
            -lboost_filesystem \
            -lboost_regex \
            -lboost_system \
            -lboost_thread \
            -lglpk

EXTRA_DIST=ARCS.py reverse_filter.py log4cxx.properties
