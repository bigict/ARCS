#ifndef _SCAFFOLD_H
#define _SCAFFOLD_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iterator>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "Edge.h"

using namespace std;

class Scaffold
{

friend ostream &operator<<(ostream &, const Scaffold &);

private:
	
	int index;
	int edge_numbers;

	vector<int > edges_of_scaffold;
	vector<int > distances_of_scaffold;

public:

	Scaffold(int, int, vector<int> &, vector<int> &);
	Scaffold();
	
	~Scaffold();
	
	//according to current edge index to get its gap distance
	int get_specific_pos_distance_of_scaffold(int) const;

	int get_edge_numbers() const;
	int get_index() const;

	vector<int> & get_edges_of_scaffold();
	vector<int> & get_distances_of_scaffold();

};

#endif
