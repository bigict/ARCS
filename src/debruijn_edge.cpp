#include "debruijn_edge.h"

#include <boost/format.hpp>

std::ostream& operator<< (std::ostream &os, const Edge &obj) {
	os << boost::format("%d\t%d\t%d\t%s") % obj.index % obj.len % obj.copy_num % obj.edge_seq;
	return os;
}
