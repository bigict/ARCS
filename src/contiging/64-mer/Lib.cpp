#include "Lib.h"
#include <string>
#include <iostream>
#include <fstream>

int Lib::get_size()
{
	return left_read_file.size();
}

void Lib::output_lib()
{
	cout << "Output libs" << endl;
	for(int i = 0; i < left_read_file.size(); i++)
	{
		cout << "LIB	" << i +1 << endl;
			cout << "	read1	" << left_read_file[i] << endl;
			cout << "	read2	" << right_read_file[i] << endl;
	}
}

void Lib::push_two_read_file(string& s1, string& s2)
{
	left_read_file.push_back(s1);
	right_read_file.push_back(s2);
}

string Lib::get_left_read_file(int i)
{
	return left_read_file[i];
}


string Lib::get_right_read_file(int i)
{
	return right_read_file[i];
}


