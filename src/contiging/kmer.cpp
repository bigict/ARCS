#include "kmer.h"

Kmer::Kmer() : _data(0) {
}

Kmer::Kmer(const std::string& seq) : _data(0) {
	sequence(seq);
}

Kmer::Kmer(const std::string& seq, size_t i, size_t j) : _data(0) {
	sequence(seq, i, j);
}

Kmer::~Kmer() {
}

const std::string Kmer::sequence() const {
	return "";
}

void Kmer::sequence(const std::string& seq, size_t i, size_t j) {
}
