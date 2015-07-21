#include "Edge.h"

Edge::Edge()
{
	//
}

Edge::~Edge()
{
	//
}

Edge::Edge(
		int _index, 
		long _len, 
		string _edge_seq, 
		int _copy_number
		)
		:index(_index),
		len(_len),
		edge_seq(_edge_seq),
		copy_number(_copy_number)

{
	nexts.resize(C_IN_EDGE,-1);
	prevs.resize(C_IN_EDGE,-1);
	//
}

ostream & operator<<( ostream &os, const Edge &obj)
{
	os << obj.index << "\t" << obj.len 	<< "\t" << obj.copy_number 	<< "\n" << obj.edge_seq ;
	return os;
}
