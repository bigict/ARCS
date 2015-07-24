#include "Split_String.h"

Split_String::Split_String(string s):str(s)
{
	//
}

vector<string> Split_String::split(const char *delim , int rep) 
{
	if (!fields.empty()) 
	{
		fields.clear();  // empty vector if necessary
	}

	string buf;
    set<char> exist;

    while(*delim != '\0'){ exist.insert(*delim);++delim;}

	for(int i = 0 ;i < str.size() ; ++i)
	{
		if (exist.find(str[i]) != exist.end())
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
