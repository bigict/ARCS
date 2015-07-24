#include "kmer.h"

Kmer::Kmer() : _data(0), _length(0) {
}

Kmer::Kmer(const std::string& seq) : _data(0), _length(0) {
    sequence(seq);
}

Kmer::Kmer(const std::string& seq, size_t i, size_t j) : _data(0), _length(0) {
    sequence(seq, i, j);
}

Kmer::~Kmer() {
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
        _length += 1;
        ++i;
    }
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

bool Kmer::operator == (const Kmer& o) const {
    return _data == o._data && _length == o._length;
}

std::ostream& operator << (std::ostream& os, const Kmer& kmer) {
    os << kmer.sequence();
    return os;
}
