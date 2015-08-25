#ifndef kmer_tbl_h_
#define kmer_tbl_h_

#include "kmer.h"
#include "component.h"

#include <iostream>
#include <tr1/unordered_map>

class Component;

class KmerTable {
public:
    typedef std::tr1::unordered_multimap< Kmer, std::pair< size_t, long >, KmerHasher> KmerList;

    KmerTable(size_t K) ;
	virtual ~KmerTable() ;

	inline size_t size() {
		return _kmerlist.size();
	}
    /*void addKmer(const Kmer& o, std::pair<size_t, size_t> pos) ;
    void filter(size_t INSERT_SIZE, const std::vector<Component>& components) ;
	std::pair<size_t, size_t> findPos(const Kmer& o) const;
    typedef std::tr1::unordered_multimap< Kmer, std::pair<size_t, size_t>, KmerHasher> KmerList;
*/
    void addKmer(const Kmer& o, std::pair<size_t, long> pos) ;
    void filter(size_t INSERT_SIZE, const ContigSet& contigset, const std::vector<Component>& components) ;
	std::pair<size_t, long> findPos(const Kmer& o) const;

	friend std::ostream& operator<<(std::ostream& os, const KmerTable& tbl) ;
private:
    KmerList _kmerlist;

private:
	size_t _k;
};

#endif
