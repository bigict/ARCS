#include "Split_String.h"


Split_String::Split_String(string s):str(s)
{
	//
}

vector<string> Split_String::split(char delim, int rep) 
{
	if (!fields.empty()) 
	{
		fields.clear();  // empty vector if necessary
	}

	string buf;
	for(int i = 0 ;i < str.size() ; ++i)
	{
		if (str[i] == delim)
		{
			if (!buf.empty())
			{
				fields.push_back(buf);
				buf.clear();
			}
		}
		else
		{
			buf += str[i];
		}
	}

	if (!buf.empty())
	{
		fields.push_back(buf);
	}
	
	return fields;
}
