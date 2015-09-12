bin_PROGRAMS=arcs
arcs_SOURCES=arff.cpp \
			arff.h \
            component.cpp \
            component.h \
            contiging.cpp \
            contiging.h \
            contigs.cpp \
            contigs.h \
            debruijn_graph.h \
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
            lp_divider.cpp \
            lp_divider.h \
            pair_read.cpp \
            pair_read.h \
            preprocess.cpp \
            preprocess.h \
            runner.h \
            remove_repeats.cpp \
            remove_repeats.h \
            scaffolding.cpp \
            scaffolding.h \
            scoring.cpp \
            scoring.h \
            uniq_edge_graph.cpp \
            uniq_edge_graph.h \
            main.cpp
arcs_CXXFLAGS=${LOG4CXX_CFLAGS}
arcs_LDFLAGS=${LOG4CXX_LIBS} \
				  -lboost_filesystem \
				  -lboost_regex \
				  -lboost_system \
				  -lglpk

EXTRA_DIST=ARCS.py log4cxx.properties
