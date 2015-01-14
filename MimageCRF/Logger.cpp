#include "Logger.h"

Log::Log(const string& filename, const bool &v){ 
	initialize(filename, v);
}

Log::Log(){}

void Log::initialize(const string& filename, const bool& v){
	verbose = v;
	logFile.open(filename.c_str());
	if(logFile.is_open() == false)
		throw("Error openning file\n");
}

void Log::initialize(const Parameters& p){
	verbose = p.verbose;
	string filename = p.outstem + ".txt";
	logFile.open(filename.c_str());
	if(logFile.is_open() == false)
		throw("Error openning file\n");
}
	
Log::~Log(){
	logFile.close();
}
