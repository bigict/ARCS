#ifndef _SPLIT_STRING_H
#define _SPLIT_STRING_H


#include<iostream>
#include<string>
#include<vector>

using namespace std;


class Split_String
{

private:

	string str;
	vector<string> fields;

public:

	Split_String(string );
	vector<string> split(char, int rep=0);

};

#endif
