#ifndef arff_h_
#define arff_h_

#include <iostream>
#include <map>
#include <string>
#include <vector>

//
// Attribute Relation File Format
//
class ARFF {
public:
    typedef std::map< std::string, std::string > Attrs;
    typedef std::vector< std::string > Item;
    typedef std::vector< Item > Data;

    Attrs attrs;
    Data data;
};

class ARFFReader {
public:
    ARFFReader(std::istream& stream) : _stream(stream) {
    }
    bool read(ARFF& arff);
private:
    std::istream& _stream;
};

#endif // arff_h_
