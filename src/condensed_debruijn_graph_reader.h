#ifndef condensed_debruijn_graph_reader_h_
#define condensed_debruijn_graph_reader_h_

#include <iostream>
#include <string>

struct CondensedDeBruijnGraphEdge {
    CondensedDeBruijnGraphEdge() : i(0), j(0), coverage(0) {
    }
    size_t i, j;
    std::string seq;
    double coverage;
};

class CondensedDeBruijnGraphReader {
public:
    CondensedDeBruijnGraphReader(std::istream& stream) : _stream(stream) {
    }
    bool read(CondensedDeBruijnGraphEdge& edge);

private:
    std::istream& _stream;
};

#endif // condensed_debruijn_graph_reader_h_
