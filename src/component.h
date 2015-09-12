#ifndef component_h_
#define component_h_

#include <iostream>
#include <string>
#include <vector>

#include "contigs.h"

class Component {
public:
    Component() : _length(0) {
    }
    virtual ~Component() {
    }

    size_t length() const {
        return _length;
    } 
    void length(size_t K, const ContigList& contigs);

    std::vector< size_t > contigs;
    std::vector< long > gaps;

private:
    friend std::ostream& operator<<(std::ostream& os, const Component& component) ;

    size_t _length;
};

typedef std::vector< Component > ComponentList;
bool ReadComponents(std::istream& stream, const ContigList& contigs, size_t K, ComponentList& components);
bool ReadComponents(const std::string& file, const ContigList& contigs, size_t K, ComponentList& components);

// component Format Specification : file name :component_iter
// Syntax
//    <component_id> :=    [0-9]+
//    <contig_id>    :=    [0-9]+ ... [0-9]+
//    <gap>          :=    [0-9]+ ... [0-9]+
// Requirements
//    1. The component id appears right after '>component' 
//    2. The number of gap must equal the number of cotigs_id minus 1

class ComponentReader {
public:
    ComponentReader(std::istream& stream) : _stream(stream) {
    }
    bool read(Component& component);

private:
    std::istream& _stream;
};

#endif // component_h_
