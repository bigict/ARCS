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
            kmer.h \
            kmer_dataset.h \
            kmer_tbl.h \
            kseq.cpp \
            kseq.h \
            preprocess.cpp \
            preprocess.h \
            runner.h \
            remove_repeats.cpp \
            remove_repeats.h \
            scaffolding.cpp \
            scaffolding.h \
            main.cpp
arcs_CXXFLAGS=${LOG4CXX_CFLAGS}
arcs_LDFLAGS=${LOG4CXX_LIBS} \
				  -lboost_filesystem \
				  -lboost_regex \
				  -lboost_system \
				  -lglpk

EXTRA_DIST=ARCS.py
