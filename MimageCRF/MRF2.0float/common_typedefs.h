/*
 * common_typedefs.h
 *
 *  Created on: May 13, 2012
 *      Author: bhole
 */

#ifndef COMMON_TYPEDEFS_H_
#define COMMON_TYPEDEFS_H_

#include <vector>
#include <string>

#ifdef __linux__
  #include <unordered_map>
  typedef std::unordered_map<std::string, int> myumap;
  typedef std::unordered_map<std::string, std::vector<float> > myumapvf;
#elif __APPLE__ && __MACH__
  #include <tr1/unordered_map>
  typedef std::tr1::unordered_map<std::string, int> myumap;
  typedef std::tr1::unordered_map<std::string, std::vector<float> > myumapvf;
#endif

typedef std::vector<float> Epair;

typedef std::vector<int> VI;
typedef std::vector<std::vector<int> > VVI;
typedef std::vector<std::vector<std::vector<int> > > VVVI;

typedef std::vector<float> VF;
typedef std::vector<std::vector<float> > VVF;
typedef std::vector<std::vector<std::vector<float> > > VVVF;
typedef std::vector<std::vector<std::vector<std::vector<float> > > > VVVVF;

typedef std::vector<double> VD;
typedef std::vector<std::vector<double> > VVD;
typedef std::vector<std::vector<std::vector<double> > > VVVD;


#endif /* COMMON_TYPEDEFS_H_ */
