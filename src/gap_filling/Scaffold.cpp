#include "Scaffold.h"

Scaffold::Scaffold(
		int _index, 
		int _edge_numbers,
		vector<int> &_edges_of_scaffold,
		vector<int> &_distances_of_scaffold
		):
	index(_index),
	edge_numbers(_edge_numbers),
	edges_of_scaffold(_edges_of_scaffold),
	distances_of_scaffold(_distances_of_scaffold)
{
	//	
}

Scaffold::Scaffold()
{

}

Scaffold::~Scaffold()
{
	//
}

int Scaffold::get_specific_pos_distance_of_scaffold( int _index) const
{
	return distances_of_scaffold[_index];
}

int Scaffold::get_edge_numbers() const
{
	return edge_numbers;
}

int Scaffold::get_index() const
{
	return index;
}

vector<int>& Scaffold::get_edges_of_scaffold()
{
	return edges_of_scaffold;
}

vector<int>& Scaffold::get_distances_of_scaffold()
{
	return distances_of_scaffold;
}



ostream &operator<<(ostream &os, const Scaffold &s)
{
	os << s.index 
		<< "\t" 
		<< s.edge_numbers 
		<< "\t";
	copy((s.edges_of_scaffold).begin(), 
			(s.edges_of_scaffold).end(),
			ostream_iterator<int>(os,","));
	os << "\t";
	copy((s.distances_of_scaffold).begin(), 
			(s.distances_of_scaffold).end(),
			ostream_iterator<int>(os,","));
	os << endl;
	
	return os;
}
