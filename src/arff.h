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
    typedef Item* ItemPtr;
    typedef std::vector< ItemPtr > Data;

    virtual ~ARFF() {
        for (Data::iterator i = data.begin(); i != data.end(); ++i) {
            delete *i;
        }
    }

    Attrs attrs;
    Data data;
};

class ARFFReader {
public:
    ARFFReader(std::istream& stream, size_t num_fields) : _stream(stream), _num_fields(num_fields) {
    }
    bool read(ARFF& arff);
private:
    std::istream& _stream;
    size_t _num_fields;
};

#endif // arff_h_
