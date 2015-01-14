/*
 * readIVMparam.cpp
 *
 *  Created on: May 22, 2010
 *      Author: bhole
 */

#include "readIVMparam.h"

using namespace boost::numeric::ublas;
using namespace std;

int load_data(const char* c, matrix<DBL_TYPE>& X, matrix<DBL_TYPE>& Y){
  ifstream dataset(c);
  if(!dataset.good()){
    return 0;
  }
  string a;
  int N = -1;
  int d = 0;
  stringstream ss (stringstream::in | stringstream::out);
  struct csv_reader* z = new csv_reader;
  ss.imbue(locale(locale(), z));
  while(!dataset.eof()){
    dataset >> a;
    ss << a;
    string t;
    d = 0;
    while(!ss.eof()){
      ss >> t;
      d++;
    }
    if(a.size() > 0)
      N = N+1;

    ss.clear();
  }

  cout << "Counting " << N << " examples ";
  cout << "of " << d-1 << " dimensions." << endl;

  dataset.clear();
  dataset.seekg(0, ios::beg);
  X = matrix<DBL_TYPE>(N,d-1);
  Y = matrix<DBL_TYPE>(N,1);
  int i = 0;
  int j;

  while(!dataset.eof()){
    dataset >> a;
    ss << a;
    string t;
    j = 0;
    if(a.size() < 1){
      break;
    }
    while(!ss.eof()){
      ss >> t;
      if(j == 0){
	(Y)(i, 0) = atof(t.c_str());
      }else{
	(X)(i, j-1) = atof(t.c_str());
      }
      j++;
    }
    i = i + 1;
    ss.clear();
    a = "";
  }
  dataset.close();
  return 1;
}

int load_one_data(const char* c, matrix<DBL_TYPE>& X){
  ifstream dataset(c);
  if(!dataset.good()){
    return 0;
  }
  string a;
  int N = -1;
  int d = 0;
  stringstream ss (stringstream::in | stringstream::out);
  struct csv_reader* z = new csv_reader;
  ss.imbue(locale(locale(), z));
  while(!dataset.eof()){
    dataset >> a;
//    cout << a << endl;
    ss << a;
    string t;
    d = 0;
    while(!ss.eof()){
      ss >> t;
//      cout << t << " ";
      d++;
    }
    if(a.size() > 0)
      N = N+1;

    ss.clear();
  }

  cout << "Counting " << N << " examples ";
  cout << "of " << d << " dimensions." << endl;

  dataset.clear();
  dataset.seekg(0, ios::beg);
  X = matrix<DBL_TYPE>(N,d);
  int i = 0;
  int j;

  while(!dataset.eof()){
    dataset >> a;
    ss << a;
    string t;
    j = 0;
    if(a.size() < 1){
      break;
    }
    while(!ss.eof()){
      ss >> t;
      (X)(i, j) = atof(t.c_str());
      j++;
    }
    i = i + 1;
    ss.clear();
    a = "";
  }
  dataset.close();
  return 1;
}

int load_one_data(const char* c, matrix<DBL_TYPE>& X, int neglect){
  ifstream dataset(c);
  if(!dataset.good()){
    return 0;
  }
  string a;
  int N = -1;
  int d = 0;
  stringstream ss (stringstream::in | stringstream::out);
  struct csv_reader* z = new csv_reader;
  ss.imbue(locale(locale(), z));
  while(!dataset.eof()){
    dataset >> a;
//    cout << a << endl;
    ss << a;
    string t;
    d = 0;
    while(!ss.eof()){
      ss >> t;
//      cout << t << " ";
      d++;
    }
    if(a.size() > 0)
      N = N+1;

    ss.clear();
  }

  if (neglect==1)
	  d = d - 1;

  cout << "Counting " << N << " examples ";
  cout << "of " << d << " dimensions." << endl;

  dataset.clear();
  dataset.seekg(0, ios::beg);
  X = matrix<DBL_TYPE>(N,d);
  int i = 0;
  int j;

  while(!dataset.eof()){
    dataset >> a;
    ss << a;
    string t;
    j = 0;
    if(a.size() < 1){
      break;
    }
    int first = 0;
    while(!ss.eof()){
      ss >> t;

      if (neglect==1 && first==0)
      {
    	  first = 1;
    	  continue;
      }

      (X)(i, j) = atof(t.c_str());
      j++;
    }
    i = i + 1;
    ss.clear();
    a = "";
  }
  dataset.close();
  return 1;
}





