bin_PROGRAMS=contiging
contiging_SOURCES= debruijn_graph.cpp \
            debruijn_graph.h \
            kmer.cpp \
            kmer.h \
            kmer_tbl.cpp \
            kmer_tbl.h \
            kseq.cpp \
            kseq.h \
            main.cpp
contiging_CXXFLAGS=${LOG4CXX_CFLAGS}
contiging_LDFLAGS=${LOG4CXX_LIBS}

SUBDIRS=32-mer 64-mer 96-mer
