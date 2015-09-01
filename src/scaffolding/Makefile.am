bin_PROGRAMS=scaffolding 
scaffolding_SOURCES=main.cpp \
            contigs.cpp	\
            contigs.h	\
            component.cpp	\
            component.h	\
            graph.cpp \
            graph.h \
            kmer.cpp	\
            kmer.h	\
            kmer_tbl.cpp	\
            kmer_tbl.h	\
            kseq.cpp \
			kseq.h	\
            pair_read.cpp	\
	        pair_read.h
scaffolding_CXXFLAGS=${LOG4CXX_CFLAGS}
scaffolding_LDFLAGS=${LOG4CXX_LIBS} \
					-lboost_filesystem \
					-lboost_regex \
				    -lboost_system
EXTRA_DIST=parse_glpk_results.py
