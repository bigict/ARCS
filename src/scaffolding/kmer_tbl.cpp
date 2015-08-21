#include "kmer_tbl.h"

KmerTable::KmerTable(size_t K) : _k(K) {
}

KmerTable::~KmerTable() {
}

void KmerTable::addKmer( const Kmer& o, std::pair<size_t, long> pos) {
    _kmerlist.insert(std::make_pair( o, pos ));
}

void KmerTable::filter(size_t INSERT_SIZE, const std::vector<Component>& components) {
    for(KmerList::iterator it=_kmerlist.begin(); it!=_kmerlist.end(); ++it) {
        size_t com_num = it->second.first;
        size_t index = it->second.second;
        size_t len = components[com_num].length();
        if (index <= 2 * INSERT_SIZE || index >= len - _k + 1 - 2 * INSERT_SIZE){
            continue;
        }
        _kmerlist.erase(it);
    }
    return ;
}
std::pair<size_t, long> KmerTable::findPos(const Kmer& o) const {
    KmerList::const_iterator it = _kmerlist.find(o);
    if (it != _kmerlist.end()) {
        return it->second;
    }
    return std::make_pair(-1, -1);
}

std::ostream& operator<<(std::ostream& os, const KmerTable& tbl) {
	os << "kmer list size:\t" << tbl._kmerlist.size() << "\t"<< "kmer list:" << std::endl;
	KmerTable::KmerList::const_iterator it = tbl._kmerlist.begin();
    while (it != tbl._kmerlist.end()) {
		os << it->first.length() << "\t" << it->first.sequence() << "\t" << it->second.first <<"\t" << it->second.second << std::endl;
		++it;
    }
	return os;
}




