#ifndef kmer_h_
#define kmer_h_

#include "utils.h"

#include <functional>
#include <string>
#include <tr1/array>
#include <tr1/unordered_map>

//
// Utilities for encodeing Nucleotide
//
class Nucleotide {
public:
    enum Code {
        Adenine  = 0x00, 
        Cytosine = 0x01, 
        Guanine  = 0x02, 
        Thymine  = 0x03,
        NUM      = 0x04
    };
    
    static int char2code(char chr) {
        chr = std::toupper(chr);
        switch (chr) {
        case 'A': return Adenine;
        case 'C': return Cytosine;
        case 'G': return Guanine;
        case 'T': return Thymine;
        }
        return Adenine;
    }
    
    static char code2char(int code) {
        switch (code) {
        case Adenine:  return 'A';
        case Cytosine: return 'C';
        case Guanine:  return 'G';
        case Thymine:  return 'T';
        }
        return 'A';
    }
};

template< size_t K >
class Kmer {
public:
    Kmer() {
    }
    Kmer(const std::string& seq) {
        sequence(seq, 0, std::string::npos);
    }
    Kmer(const std::string& seq, size_t i, size_t j) {
        sequence(seq, i, j);
    }
    Kmer(const Kmer< K >& o) {
        _data = o._data;
    }

    static size_t length() {
        return _K;
    }
    static void length(size_t L) {
        _K = L;
    }
    void sequence(const std::string& seq, size_t i, size_t j) {
        j = (j == std::string::npos || j >= seq.length()) ? seq.length() : j;

        std::fill(_data.begin(), _data.end(), 0);

        size_t k = 0;
        while (i < j && k < Kmer< K >::_K) {
            size_t code = Nucleotide::char2code(seq[i]);
            size_t m = (2 * k) / SIZEOF_BITS(size_t), n = SIZEOF_BITS(size_t) - ((2 * k) % SIZEOF_BITS(size_t)) - 2;
            _data[m] |= (code << n);
            ++i, ++k;
        }
    }
    std::string sequence() const {
        std::string seq;

        size_t k = 0;
        while (k < Kmer< K >::_K) {
            size_t m = (2 * k) / SIZEOF_BITS(size_t), n = SIZEOF_BITS(size_t) - ((2 * k) % SIZEOF_BITS(size_t)) - 2;
            size_t code = (_data[m] >> n) & 0x03;
            seq += Nucleotide::code2char(code);
            ++k;
        }

        return seq;
    }
    char nucleotide(size_t k) const {
        size_t m = (2 * k) / SIZEOF_BITS(size_t), n = SIZEOF_BITS(size_t) - ((2 * k) % SIZEOF_BITS(size_t)) - 2;
        size_t code = (_data[m] >> n) & 0x03;
        return Nucleotide::code2char(code);
    }

    void shift(char nucleotide) {
        shift(Nucleotide::char2code(nucleotide));
    }
    void shift(int nucleotide) {
        for (size_t i = 1; i < _data.size(); ++i) {
            _data[i - 1] <<= 2;
            _data[i - 1] |= (_data[i] >> (SIZEOF_BITS(size_t) - 2) & 0x03);
        }
        if (!_data.empty()) {
            size_t k = Kmer< K >::_K - 1;
            size_t m = (2 * k) / SIZEOF_BITS(size_t), n = SIZEOF_BITS(size_t) - ((2 * k) % SIZEOF_BITS(size_t)) - 2;
            _data[m] <<= 2;
            _data[m] |= (((size_t)nucleotide) << n);
        }
    }

    void unshift(char nucleotide) {
        unshift(Nucleotide::char2code(nucleotide));
    }
    void unshift(int nucleotide) {
        if (!_data.empty()) {
            size_t k = Kmer< K >::_K - 1;
            size_t m = (2 * k) / SIZEOF_BITS(size_t), n = SIZEOF_BITS(size_t) - ((2 * k) % SIZEOF_BITS(size_t)) - 2;
            _data[m] >>= n + 2;
            _data[m] <<= n + 2;
        }
        for (size_t i = _data.size(); i > 1; --i) {
            _data[i - 1] >>= 2;
            _data[i - 1] |= ((_data[i - 2] & 0x03) << (SIZEOF_BITS(size_t) - 2));
        }
        if (!_data.empty()) {
            _data[0] >>= 2;
            _data[0] |= (((size_t)nucleotide) << (SIZEOF_BITS(size_t) - 2));
        }
    }

    Kmer< K > operator + (const char c) const {
        Kmer kmer = *this;
        kmer.shift(c);
        return kmer;
    }

    Kmer< K > operator - (const char c) const {
        Kmer kmer = *this;
        kmer.unshift(c);
        return kmer;
    }

    bool operator < (const Kmer< K >& o) const {
        return _data < o._data;
    }
    bool operator > (const Kmer< K >& o) const {
        return _data > o._data;
    }
    bool operator == (const Kmer< K >& o) const {
        return _data == o._data;
    }
    bool operator != (const Kmer< K >& o) const {
        return _data != o._data;
    }

    size_t hash() const {
        static std::hash< size_t > hasher;
        return hasher(_data[0]);
    }

private:
    std::tr1::array< size_t,  (2 * K + SIZEOF_BITS(size_t) - 1) / SIZEOF_BITS(size_t) > _data;

    static size_t _K;
};

template< size_t K >
size_t Kmer< K >::_K = 31;

template< size_t K >
struct KmerHasher {
    size_t operator()(const Kmer< K >& o) const {
        return o.hash();
    }
};

template< size_t K >
std::ostream& operator << (std::ostream& os, const Kmer< K >& kmer) {
    os << kmer.sequence();
    return os;
}

template< size_t K, class V >
using KmerTable = std::tr1::unordered_map< Kmer< K >, V, KmerHasher< K > >;

#endif // kmer_h_
