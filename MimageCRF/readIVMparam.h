/*
 * readIVMparam.h
 *
 *  Created on: May 22, 2010
 *      Author: bhole
 */

#ifndef READIVMPARAM_H_
#define READIVMPARAM_H_


#include <iterator>
#include <string>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/assign/std/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <boost/timer.hpp>
#include "ivm.hpp"
//#include "klropts.hpp"
//#include "cblas_mm.hpp"



struct csv_reader: std::ctype<char> {
  csv_reader(): std::ctype<char>(get_table()) {}
  static std::ctype_base::mask const* get_table() {
    static std::vector<std::ctype_base::mask> rc(table_size, std::ctype_base::mask());

    rc[','] = std::ctype_base::space;
    rc['\n'] = std::ctype_base::space;
    return &rc[0];
  }
};

struct tab_reader: std::ctype<char> {
  tab_reader(): std::ctype<char>(get_table()) {}
  static std::ctype_base::mask const* get_table() {
    static std::vector<std::ctype_base::mask> rc(table_size, std::ctype_base::mask());

    rc['\t'] = std::ctype_base::space;
    rc['\n'] = std::ctype_base::space;
    return &rc[0];
  }
};

using namespace boost::numeric::ublas;

int load_data(const char* c, matrix<DBL_TYPE>& X, matrix<DBL_TYPE>& Y);
int load_one_data(const char* c, matrix<DBL_TYPE>& X);
int load_one_data(const char* c, matrix<DBL_TYPE>& X, int neglect);

#endif /* READIVMPARAM_H_ */
