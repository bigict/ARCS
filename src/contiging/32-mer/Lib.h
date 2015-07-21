#ifndef _lib_
#define _lib_
#include <vector>
#include <string>

using namespace std;

class Lib
{
public:
	int get_size();
	void push_two_read_file(string&, string&);
	
	string get_left_read_file(int);
	string get_right_read_file(int);
	
	void output_lib();
	
private:
	vector<string> left_read_file;
	vector<string> right_read_file;
};

#endif
