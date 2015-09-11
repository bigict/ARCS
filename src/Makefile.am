bin_PROGRAMS=arcs
arcs_SOURCES=arff.cpp \
			arff.h \
            contiging.cpp \
            contiging.h \
            debruijn_graph.h \
            kmer.h \
            kmer_dataset.h \
            kmer_tbl.h \
            kseq.cpp \
            kseq.h \
            preprocess.cpp \
            preprocess.h \
            runner.h \
            scaffolding.cpp \
            scaffolding.h \
            main.cpp
arcs_CXXFLAGS=${LOG4CXX_CFLAGS}
arcs_LDFLAGS=${LOG4CXX_LIBS} \
				  -lboost_filesystem \
				  -lboost_regex \
				  -lboost_system

#SUBDIRS=contiging scaffolding divideLP remove_repeats gap_filling
EXTRA_DIST=ARCS.py
