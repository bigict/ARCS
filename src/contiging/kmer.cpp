#include "kmer.h"

#include <boost/foreach.hpp>

Kmer::Kmer() : _data(0), _length(0) {
}

Kmer::Kmer(const std::string& seq) : _data(0), _length(0) {
    sequence(seq);
}

Kmer::Kmer(const Kmer& o) {
    _data = o._data;
    _length = o._length;
}

Kmer::Kmer(const std::string& seq, size_t i, size_t j) : _data(0), _length(0) {
    sequence(seq, i, j);
}

Kmer::~Kmer() {
}

Kmer Kmer::subKmer(size_t i, size_t j) const {
    j = (j == -1) ? _length : j;
    bigint musk(0x03); 
    
    Kmer kmer;
    while (i < j) {
        int code = boost::lexical_cast< int >((_data >> ((_length - i - 1) * 2)) & musk);
        kmer._data = (kmer._data << 2) + code;
        kmer._length += 1;
        ++i;
    }
    return kmer;
}

Nucleotide::Code Kmer::nucleotide(size_t i) const {
    BOOST_ASSERT(i >=0 && i < _length);
    if (i >= 0 && i < _length) {
        bigint musk(0x03); 
        return (Nucleotide::Code)boost::lexical_cast< int >((_data >> ((_length - i - 1) * 2)) & musk);
    }
    return Nucleotide::Adenine;
}

const std::string Kmer::sequence() const {
    std::string seq;

    bigint musk(0x03); 
    for (size_t i = 0; i < _length; ++i) {
        int code = boost::lexical_cast< int >((_data >> ((_length - i - 1) * 2)) & musk);
        seq += Nucleotide::code2char(code);
    }

    return seq;
}

void Kmer::sequence(const std::string& seq, size_t i, size_t j) {
    j = (j == std::string::npos || j >= seq.length()) ? seq.length() : j;

    _data = 0;
    _length = 0;

    while (i < j) {
        int code = Nucleotide::char2code(seq[i]);
        _data  = (_data << 2) + code;
        ++_length;
        ++i;
    }
}

Nucleotide::Code Kmer::pop() {
    if (_length > 0) {
        bigint musk(0x03); 
        int code = boost::lexical_cast< int >((_data >> ((_length - 1) * 2)) & musk);
        --_length;
        return (Nucleotide::Code)code;
    }
    return Nucleotide::Adenine;
}

void Kmer::push(Nucleotide::Code code) {
    _data  = (_data << 2) + (int)code;
    ++_length;
}

Kmer& Kmer::operator = (const std::string& seq) {
    sequence(seq);
    return *this;
}

Kmer& Kmer::operator = (const Kmer& o) {
    if (&o != this) {
        _data = o._data;
        _length = o._length;
    }
    return *this;
}

Kmer Kmer::operator + (const char c) {
    Kmer kmer = *this;
    kmer += c;
    return kmer;
}

Kmer Kmer::operator + (const std::string& seq) {
    Kmer kmer = *this;

    BOOST_FOREACH(const char c, seq) {
        kmer += c;
    }

    return kmer;
}

Kmer Kmer::operator + (const Kmer& o) {
    Kmer kmer = *this;

    kmer._data <<= o._length;
    kmer._data += o._data;
    kmer._length += o._length;

    return kmer;
}

Kmer& Kmer::operator += (const char c) {
    int code = Nucleotide::char2code(c);
    _data  = (_data << 2) + code;
    ++_length;

    return *this;
}

Kmer& Kmer::operator += (const std::string& seq) {
    BOOST_FOREACH(const char c, seq) {
       *this += c;
    }
    return *this;
}

Kmer& Kmer::operator += (const Kmer& o) {
    if (o._length > 0) {
            _data <<= o._length;
            _data += o._data;
            _length += o._length;
    }
    return *this;
}


bool Kmer::operator == (const Kmer& o) const {
    return _data == o._data && _length == o._length;
}

std::ostream& operator << (std::ostream& os, const Kmer& kmer) {
    os << kmer.sequence();
    return os;
}

