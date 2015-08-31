#ifndef kmer_tbl_h_
#define kmer_tbl_h_

#include "component.h"
#include "contigs.h"
#include "kmer.h"

#include <iostream>
#include <tr1/unordered_map>

class Component;

class KmerTable {
public:
    KmerTable(size_t K) : _K(K) {
    }
	virtual ~KmerTable() {
    }

	size_t size() const {
		return _kmerlist.size();
	}

    void addKmer(const Kmer& o, std::pair<size_t, long> pos) ;
	std::pair<size_t, long> findPos(const Kmer& o) const;

private:
	friend std::ostream& operator<<(std::ostream& os, const KmerTable& tbl) ;

    typedef std::tr1::unordered_multimap< Kmer, std::pair<size_t, long>, KmerHasher> KmerList;
    KmerList _kmerlist;

private:
	size_t _K;
};

typedef std::pair< size_t, long > KmerPosition;
typedef std::tr1::unordered_multimap< Kmer, KmerPosition, KmerHasher > KmerList;

size_t BuildKmerTable(size_t K, size_t insert_size, const ContigList& contigs, const ComponentList& components, KmerList& tbl);

#endif // kmer_tbl_h_

