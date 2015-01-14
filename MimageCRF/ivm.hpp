/*
 * ivm.hpp
 *
 *  Created on: May 22, 2010
 *      Author: bhole
 */

#ifndef IVM_HPP_
#define IVM_HPP_


#include <complex>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <math.h>
//#include "cblas_mm.hpp"

#include "common.h"

#define INF std::numeric_limits<DBL_TYPE>::max();
#define temp_type_v ublas::vector<DBL_TYPE>
#define EDRPS 1e-8
#define disp(x) (std::cout << x << std::endl)
#define DBL_TYPE double

//enum kernelType { RBF, CHI2, POLY, LINEAR, EUCLIDEAN, NONE };

using namespace boost::numeric;

class IVM {

  ublas::matrix<DBL_TYPE>* X;
  ublas::matrix<DBL_TYPE>* w;
  kernelType kernel;
  DBL_TYPE lambda;
  DBL_TYPE p1;

public:

  int K;

  ~IVM(){
  }

  ublas::matrix<DBL_TYPE> classify(ublas::matrix<DBL_TYPE>* Xtest){

	int Ntr = X->size1();
    int N = Xtest->size1();

//    std::cout << " inff = " << Ntr << "  " << N << "  " << X->size2() <<std::endl;


	DBL_TYPE xi = 1;
	ublas::matrix<DBL_TYPE> f(N, K);


	for(int l = 0; l < N; l++){
		ublas::matrix_row< ublas::matrix<DBL_TYPE> > v1(*Xtest, l);

		ublas::vector< DBL_TYPE > u(v1.size());

		f(l, K-1) = 0;

		for (int mm=0; mm<K-1; mm++) {

			f(l, mm) = 0;

			for(int ltr = 0; ltr < Ntr; ltr++){
				ublas::matrix_row< ublas::matrix<DBL_TYPE> > v2(*X, ltr);

				switch(kernel){
					case RBF:
						u = v1 - v2;

						f(l, mm) += (*w)(ltr, mm) * (1/xi)*exp(-inner_prod(u, u)/(2.0*pow(p1,2)));
						break;
					case POLY:
						f(l, mm) += (*w)(ltr, mm) * (1/xi) * (pow(inner_prod(v1, v2) + 1, p1));
						break;
					case LINEAR:
					default:
						f(l, mm) += inner_prod(v1, v2);
				}
			}
		}
	}

    return f;
  }


  ublas::vector<DBL_TYPE> classify_one(ublas::vector<DBL_TYPE>* Xtest, int pass, int pass2){

	int Ntr = X->size1();
    int N = Xtest->size();

 //   std::cout << " inff = " << Ntr << "  " << N << "  " << X->size2() <<std::endl;

	DBL_TYPE xi = 1;
	ublas::vector<DBL_TYPE> f(K);

	ublas::vector< DBL_TYPE > u(Xtest->size());

	f(K-1) = 0;

	for (int mm=0; mm<K-1; mm++) {

		f(mm) = 0;

		for(int ltr = 0; ltr < Ntr; ltr++){
			ublas::matrix_row< ublas::matrix<DBL_TYPE> > v2(*X, ltr);

			switch(kernel){
				case RBF:

					if (pass==1) {
						for (int jki=0; jki<N; jki++)
							std::cout << " " << (*Xtest)(jki);
						std::cout << std::endl;
						for (int jki=0; jki<N; jki++)
							std::cout << " " << v2(jki);
						std::cout << std::endl;

					}


					u = *Xtest - v2;

					f(mm) += (*w)(ltr, mm) * (1/xi)*exp(-inner_prod(u, u)/(2.0*pow(p1,2)));

					if (pass==1)
						std::cout<< " mm = " << mm << " ltr = " << ltr  << " f val = " << f(mm) << std::endl;

					break;
				case POLY:
					f(mm) += (*w)(ltr, mm) * (1/xi) * (pow(inner_prod(*Xtest, v2) + 1, p1));
					break;
				case LINEAR:
				default:
					f(mm) += inner_prod(*Xtest, v2);
			}
		}
	}


    return f;
  }


  DBL_TYPE classify_one(ublas::vector<DBL_TYPE>* Xtest, int j) {

 	int N = Xtest->size();

 	DBL_TYPE xi = 1;
 	DBL_TYPE Kval;

 	ublas::vector< DBL_TYPE > u(Xtest->size());

	ublas::matrix_row< ublas::matrix<DBL_TYPE> > v2(*X, j);

	switch(kernel){
		case RBF:
			u = *Xtest - v2;

			Kval = (1/xi)*exp(-inner_prod(u, u)/(2.0*pow(p1,2)));
			break;
		case POLY:
			Kval = (1/xi) * (pow(inner_prod(*Xtest, v2) + 1, p1));
			break;
		case LINEAR:
		default:
			Kval = inner_prod(*Xtest, v2);
	}


     return Kval;
   }



  ublas::vector<DBL_TYPE> classify_klrexp(ublas::vector<DBL_TYPE>* Xtest){

	int Ntr = X->size1();
    int N = Xtest->size(); // dimension of relevance vector

	DBL_TYPE xi = 1;
	ublas::vector<DBL_TYPE> f(K);

	ublas::vector< DBL_TYPE > u(Xtest->size());

	for (int mm=0; mm<K; mm++) {

/*		f(mm) = +0.448545; // -rho for first class
		if (mm==1)
			f(mm) = -f(mm);
*/
		f(mm) = 0;

		for(int ltr = 0; ltr < Ntr; ltr++){
			ublas::matrix_row< ublas::matrix<DBL_TYPE> > v2(*X, ltr);

			switch(kernel){
				case RBF: {
					u = *Xtest - v2;
					f(mm) += (*w)(ltr, mm) * (1/xi)*exp(-inner_prod(u, u)/(2.0*pow(p1,2)));
					break;
				}
				case POLY:
					f(mm) += (*w)(ltr, mm) * (1/xi) * (pow(inner_prod(*Xtest, v2) + 1, p1));
					break;
				case LINEAR:
				default:
					f(mm) += inner_prod(*Xtest, v2);
			}
		}
//		std::cout << "\n";
		// temporary statement for svm
		//f(mm) = -f(mm);
//		std::cout << f(mm) << " ";
	}
//	std::cout << "\n";

    return f;
  }


  // this return K(X_t, x_j)
  // remember no w multiplication because derivative removes it
  DBL_TYPE classify_klrexp(ublas::vector<DBL_TYPE>* Xtest, int j) {

 	int N = Xtest->size();

 	DBL_TYPE xi = 1;
 	DBL_TYPE Kval;

 	ublas::vector< DBL_TYPE > u(Xtest->size());

	ublas::matrix_row< ublas::matrix<DBL_TYPE> > v2(*X, j);

	//std::cout<<"  v2 "<<v2<<std::endl;

	switch(kernel){
		case RBF: {
			u = *Xtest - v2;
			//std::cout<<" u  "<<u<<std::endl;

			Kval = (1/xi)*exp(-inner_prod(u, u)/(2.0*pow(p1,2)));
			break;
		}
		case POLY:
			Kval = (1/xi) * (pow(inner_prod(*Xtest, v2) + 1, p1));
			break;
		case LINEAR:
		default:
			Kval = inner_prod(*Xtest, v2);
	}


     return Kval;
   }




  IVM(ublas::matrix<DBL_TYPE>* Xt,
	  ublas::matrix<DBL_TYPE>* wparam,
      kernelType Kernel,
      DBL_TYPE Lambda,
      int k,
      DBL_TYPE P) {

	X = Xt;
	w = wparam;

    kernel = Kernel;
    lambda = Lambda;
    p1 = P;
    K = k; // number of classes
  }

};


#endif /* IVM_HPP_ */
