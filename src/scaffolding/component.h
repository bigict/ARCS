#ifndef component_h_
#define component_h_

#include <string>
#include <iostream>

#include "contigs.h"
#include "kmer_tbl.h"

class KmerTable; 

class Component {
public:
    Component(size_t K) : _k(K) {
    }
    virtual ~Component() {
    }
//  size_t produceKmer(const ContigSet& c, KmerTable& tbl, size_t component_no);
    size_t produceKmerForInsertSize(const ContigSet& c, KmerTable& tbl, size_t component_no); 
    size_t produceKmerForPairRead(const ContigSet& c, KmerTable& tbl, size_t component_no, size_t INSERT_SIZE);

    size_t getLen() const { return _len; } 
    void initializeLen(const ContigSet& c);
    void reset() {
        _contig_id.clear();
        _gap.clear();
    }

    std::vector< size_t > _contig_id;
    std::vector< long > _gap;

private:
    friend std::ostream& operator<<(std::ostream& os, const Component& component) ;

    size_t _len;
    size_t _k;
    friend class ComponentReader;
};

class ComponentReader {
public:
    ComponentReader(std::istream& stream) : _stream(stream) {
    }

    bool read(Component& component);

private:
    std::istream& _stream;
};

#endif // component_h_
