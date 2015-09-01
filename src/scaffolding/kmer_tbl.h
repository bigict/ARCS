#ifndef kmer_tbl_h_
#define kmer_tbl_h_

#include "component.h"
#include "contigs.h"
#include "kmer.h"

#include <iostream>
#include <tr1/unordered_map>

typedef std::pair< size_t, long > KmerPosition;
typedef std::tr1::unordered_multimap< Kmer, KmerPosition, KmerHasher > KmerList;

size_t BuildKmerTable_insertsize(size_t K, size_t insert_size, const ContigList& contigs, const ComponentList& components, KmerList& tbl);
size_t BuildKmerTable_pairends(size_t K, size_t insert_size, size_t edge_cutoff, const ContigList& contigs, const ComponentList& components, KmerList& tbl);

#endif // kmer_tbl_h_

