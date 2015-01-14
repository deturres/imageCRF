#ifndef LOGGER_FILE
#define LOGGER_FILE

#include <fstream>
#include <string>
#include <iostream>
#include "Parameters.h"
#include <vector>

using std::cout;
using std::ofstream;
using std::string;
using std::vector;

class Log{
public:
	Log();
	Log(const string&, const bool &v);
	void initialize(const string&, const bool&);
	void initialize(const Parameters&);

	~Log();

private:
	ofstream logFile;
	bool verbose;

public:


	// compiler does not like template function outside of class.
	template <class t>
	Log& operator<<(const vector<t>& a){
		if(verbose == true){
			for(vector<t>::const_iterator b = a.begin(); b != a.end(); ++b)
				*this << *b << " ";
		}
		return *this;
	}

	// compiler does not like template function outside of class.
	template <class t>
	Log& operator<<(const t& a){
		if(verbose == true){
			cout << a;
			logFile << a;
		}
		return *this;
	}
};


#endif