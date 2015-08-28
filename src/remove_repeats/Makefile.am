bin_PROGRAMS=remove_repeats 
remove_repeats_SOURCES=main.cpp \
            component.cpp \
            component.h \
            uniq_edge_graph.cpp \
            uniq_edge_graph.h
remove_repeats_CXXFLAGS=${LOG4CXX_CFLAGS}
remove_repeats_LDFLAGS=${LOG4CXX_LIBS} \
					   -lboost_filesystem \
					   -lboost_regex \
					   -lboost_system
