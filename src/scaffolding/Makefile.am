bin_PROGRAMS=scaffolding
scaffolding_SOURCES=main.cpp \
            linearization.cpp	\
            linearization.h	\
            pos_recorder.cpp	\
            pos_recorder.h	\
            uniq_edge_graph.cpp	\
            uniq_edge_graph.h	\
            kmer.h	\
            kmer.cpp	\
            kmer_hash.h	\
            kmer_hash.cpp
scaffolding_CXXFLAGS=${LOG4CXX_CFLAGS}
scaffolding_LDFLAGS=${LOG4CXX_LIBS} \
            -lgmp
LIBS += -lpthread
