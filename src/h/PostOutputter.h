#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

class PostOutputter
{
private:

//!	File stream for output
	ofstream OutputFile;

protected:

//!	Constructor
	PostOutputter(string FileName);

//!	Designed as a single instance class
	static PostOutputter* _instance;

public:

//!	Return pointer to the output file stream
	ofstream* GetOutputFile() { return &OutputFile; }

//!	Return the single instance of the class
	static PostOutputter* Instance(string FileName = " ");


//!	Output element stresses 
	void OutputElementStress();

//! Overload the operator <<
	template <typename T>
	PostOutputter& operator<<(const T& item) 
	{
		OutputFile << item;
		return *this;
	}

	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	PostOutputter& operator<<(CharOstream& (*op)(CharOstream&)) 
	{
		op(OutputFile);
		return *this;
	}
};