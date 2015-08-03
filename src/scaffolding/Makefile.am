bin_PROGRAMS=scaffolding
scaffolding_SOURCES=main.cpp \
            Edge.h	\
            Linearization.cpp	\
            Linearization.h	\
            Pos_Recorder.cpp	\
            Pos_Recorder.h	\
            Uniq_Edge_Graph.cpp	\
            Uniq_Edge_Graph.h	\
            Kmer.h	\
            Kmer.cpp	\
            Kmer_Hash.h	\
            Kmer_Hash.cpp \
            Lib.h \
            Lib.cpp 
scaffolding_CXXFLAGS=${LOG4CXX_CFLAGS}
scaffolding_LDFLAGS=${LOG4CXX_LIBS} \
            -lgmp
LIBS += -lpthread
