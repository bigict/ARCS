#ifndef _FUNC_H
#define _FUNC_H

#include <iostream>
#include <string>

using namespace std;

inline string int2Str(long number )
{
	char temp[100];
	bool flag = false;
	if(number < 0)
	{
		number = -number;
		flag = true;
	}

	int digits = 0;
	int t_number = number;
	do
	{
		digits++;
		t_number /= 10;
	}while(t_number);
	
	int index = 0;
	if(flag)
	{
		temp[index++] = '-';
	}
	
	for(int i = digits+index-1; i >= index ; --i )
	{
		temp[i] = (char)(number%10 + '0');
		number/=10;
	}
	temp[digits+index] = '\0';
	return string(temp);
}

#endif
