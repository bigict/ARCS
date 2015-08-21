#ifndef component_h_
#define component_h_

#include <string>
#include <iostream>

#include "contigs.h"
#include "kmer_tbl.h"

class KmerTable; 

class Component {
public:
    Component(size_t K) : _K(K) {
    }
    virtual ~Component() {
    }
//  size_t produceKmer(const ContigSet& c, KmerTable& tbl, size_t component_no);
    size_t produceKmerForInsertSize(const ContigSet& c, KmerTable& tbl, size_t component_no); 
    size_t produceKmerForPairRead(const ContigSet& c, KmerTable& tbl, size_t component_no, size_t INSERT_SIZE);

    size_t length() const { return _length; } 
    void initializeLen(const ContigSet& c);
    void reset() {
        _contig_id.clear();
        _gap.clear();
    }

    std::vector< size_t > _contig_id;
    std::vector< long > _gap;

    static size_t length(size_t K, const ContigSet& contigs, const Component& component);
private:
    friend std::ostream& operator<<(std::ostream& os, const Component& component) ;
    friend class ComponentReader;

    size_t _length;
    size_t _K;
};

class ComponentReader {
public:
    ComponentReader(std::istream& stream) : _stream(stream) {
    }

    bool read(Component& component);

private:
    std::istream& _stream;
};

typedef std::vector< Component > ComponentList;
bool ReadComponents(std::istream& stream, ComponentList& components);
bool ReadComponents(const std::string& filename, ComponentList& components);

#endif // component_h_
