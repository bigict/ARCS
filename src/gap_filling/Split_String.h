#ifndef _SPLIT_STRING_H
#define _SPLIT_STRING_H

#include <bits/stdc++.h>

using namespace std;


class Split_String
{

private:

	string str;
	vector<string> fields;

public:

	Split_String(string );
	vector<string> split(const char *, int rep=0);
};

#endif
