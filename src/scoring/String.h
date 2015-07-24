#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>

using namespace std;

class String
{
	friend ostream &operator<<(ostream &os, const String &);
public:
	String(size_t size = 0);
	String(const char* str);
	String(const String &str);
	String& operator=(const String &str);
	String& operator+=(const String &str);
	String operator+(const String &str);
	char &operator[](size_t index){return data[index];}
	char operator[](size_t index) const{return data[index];}
	
	char &at(size_t index);
	char at(size_t index) const;

	const char *c_str() const{return data;}
	~String();

	size_t size() const{return len;}
	bool empty() const {return (len == 0);}
	void clear() {len = 0; data[0] == '\0';}
private:
	char *data;
	size_t len;
	size_t capacity;

};

String::String(size_t size)
{
	len = size;
	data = new char[len+1];
	data[0] = '\0';
	capacity = len + 1;
}

String::String(const char *str)
{
	len = strlen(str);
	data = new char[len+1];
	strcpy(data, str);
	capacity = len+1;
}
String::String(const String &str)
{
	data = new char[str.len+1];
	strcpy(data,str.data);
	len = str.len;
	capacity = len+1;
}
String::~String()
{
	delete [] data;
}
String & String::operator=(const String &str)
{
	if (str.data != data)
	{
		if(data)
		{
			delete [] data;
		}
		data = new char[str.len+1];
		strcpy(data,str.data);
		len = str.len;
		capacity = len + 1;
	}
	return *this;
}
String & String::operator+=(const String &str)
{
	if(capacity >= len+str.len+1)
	{
		strcpy(data,str.data);
		len = len + str.len;
		return *this;
	}

	char *sum = new char[len+str.len+1];
	strcpy(sum,data);
	strcpy(sum+len,str.data);
	delete [] data;
	data = sum;
	len = len+str.len;
	capacity = len+1;
	return *this;
}

String String::operator+(const String &str)
{
	String sum = *this;
	sum += str;
	return sum;
}

char & String::at(size_t index)
{
	if(index <0 || index > len)
	{
		throw runtime_error("out of range in String");
		exit(1);
	}
	return data[index];
}

char String::at(size_t index) const
{
	if(index <0 || index > len)
	{
		throw runtime_error("out of range in String");
		exit(1);
	}
	return data[index];

}
ostream &operator<<(ostream &os, const String &obj)
{
	os << obj.data;
	return os;
}


