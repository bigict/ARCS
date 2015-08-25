#ifndef component_h_
#define component_h_

#include <string>
#include <iostream>

#include "contigs.h"
#include "kmer_tbl.h"

class KmerTable; 

class Component {
public:
    Component() {
    }
    virtual ~Component() {
    }
//  size_t produceKmer(const ContigSet& c, KmerTable& tbl, size_t component_no);
    size_t produceKmerForInsertSize(size_t K, const ContigSet& c, KmerTable& tbl, size_t component_no); 
    size_t produceKmerForPairRead(size_t K, const ContigSet& c, KmerTable& tbl, size_t component_no, size_t INSERT_SIZE);

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
