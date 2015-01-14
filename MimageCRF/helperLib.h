#ifndef HELPERLIB_H
#define HELPERLIB_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include "crfmodel.h"
#include "expectations.h"
#include "evaldisps.h"


#ifdef _MSC_VER        // Visual C++
#include "dirent.h"
const char dirs = '\\';
#else
const char dirs = '/';
# include <dirent.h>
#endif

using std::vector;

// read images from a directory.
int readImages(char *, std::vector<CByteImage>&, bool);
int readImages(char *, std::vector<CImage2>&, bool);
int readDirectories(std::string& , std::vector<string>&, std::vector<string>&, bool testDir );
int readImages(string, std::vector<CImage2>&);
int readImages(string, std::vector<CByteImage>&);
int readImages(string, int, std::vector<CImage2>&);

int readImagesPaths(string dirName, std::vector<string>& filepaths);

#endif
