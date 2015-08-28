#ifndef component_h__
#define component_h__

#include <string>
#include <iostream>


#include "contigs.h"
#include "kmer_tbl.h"

class KmerTable; 

/*struct element_component {
     size_t contig_id;
     int gap; 
};*/

class Component {
public:
    explicit Component(size_t K) : _k(K) {
    }

    virtual ~Component() {}
    //size_t produceKmer(const ContigSet& c, KmerTable& tbl, size_t component_no);
    size_t produceKmerForInsertSize(const ContigSet& c, KmerTable& tbl, size_t component_no); 
    size_t produceKmerForPairRead(const ContigSet& c, KmerTable& tbl, size_t component_no, size_t INSERT_SIZE);
    inline size_t getLen() const { return _len; } 
    void initializeLen(const ContigSet& c);
    friend std::ostream& operator<<(std::ostream& os, const Component& component) ;
    
    std::vector< size_t > _contig_id;
    std::vector< long > _gap;
    //std::vector< element_component > _component;

private:
    size_t _len;
    size_t _k;
    friend class ComponentReader;
};

class ComponentReader {
public:
    ComponentReader(std::istream& stream);
    bool read(Component& component);

private:
    std::istream& _stream;
};

// component Format Specification : file name :component_iter
// Syntax
//    <component_id> :=    [0-9]+
//    <contig_id>    :=    [0-9]+ ... [0-9]+
//    <gap>          :=    [0-9]+ ... [0-9]+
// Requirements
//    1. The component id appears right after '>component' 
//    2. The number of gap must equal the number of cotigs_id minus 1


#endif // component_h_
