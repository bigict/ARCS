#ifndef component_h_
#define component_h_

#include <iostream>
#include <string>
#include <vector>

struct Component {
    struct Item {
        Item(size_t contig, int gap) : contig(contig), gap(gap) {
        }
        size_t contig;
        int gap;
    };
    typedef std::vector< size_t > ContigList;
    typedef std::vector< int > GapList;

    void init(const ContigList& contigs, const GapList& gaps);

    typedef std::vector< Item > ItemList;
    ItemList items;
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
