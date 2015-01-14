/*
 * cblas_mm.hpp
 *
 *  Created on: May 22, 2010
 *      Author: bhole
 */

#ifndef CBLAS_MM_HPP_
#define CBLAS_MM_HPP_

#include <iostream>
#include <iomanip>
extern "C"{
#include "cblas.h"
}

#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>
#include <boost/numeric/ublas/detail/matrix_assign.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>


using namespace boost::numeric;
using namespace std;

BOOST_UBLAS_INLINE
void
cblas_prod (ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& C,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& A,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& B) {
  // C(m,n) <- aplha A(m,k) B(k,n) + beta C(m,n)
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  int m = C.size1();
  int n = C.size2();
  int k = A.size2();
  BOOST_UBLAS_CHECK (C.size1() == A.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (A.size2() == B.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (B.size2() == C.size2(), ublas::bad_size ());
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	      m, n, k, 1.0, &A.data()[0], A.size2(), &B.data()[0], B.size2(), 1.0, &C.data()[0], C.size2());
}

BOOST_UBLAS_INLINE
void
cblas_prod (ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& C,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& A,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& B) {
  // C(m,n) <- aplha A(m,k) B(k,n) + beta C(m,n)
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  int m = C.size1();
  int n = C.size2();
  int k = A.size2();
  BOOST_UBLAS_CHECK (C.size1() == A.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (A.size2() == B.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (B.size2() == C.size2(), ublas::bad_size ());
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
	      m, n, k, 1.0, &A.data()[0], A.size1(), &B.data()[0], B.size2(), 1.0, &C.data()[0], C.size2());
}

BOOST_UBLAS_INLINE
void
cblas_prod (ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& C,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& A,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& B) {
  // C(m,n) <- aplha A(m,k) B(k,n) + beta C(m,n)
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  int m = C.size1();
  int n = C.size2();
  int k = A.size2();
  BOOST_UBLAS_CHECK (C.size1() == A.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (A.size2() == B.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (B.size2() == C.size2(), ublas::bad_size ());
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	      m, n, k, 1.0, &A.data()[0], A.size2(), &B.data()[0], B.size1(), 1.0, &C.data()[0], C.size2());
}

BOOST_UBLAS_INLINE
void
cblas_prod (ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& C,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& A,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& B) {
  // C(m,n) <- aplha A(m,k) B(k,n) + beta C(m,n)
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  int m = C.size1();
  int n = C.size2();
  int k = A.size2();
  BOOST_UBLAS_CHECK (C.size1() == A.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (A.size2() == B.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (B.size2() == C.size2(), ublas::bad_size ());
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
	      m, n, k, 1.0, &A.data()[0], A.size1(), &B.data()[0], B.size1(), 1.0, &C.data()[0], C.size2());
}

BOOST_UBLAS_INLINE
void
cblas_prod (ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& C,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& A,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& B) {
  // C(m,n) <- aplha A(m,k) B(k,n) + beta C(m,n)
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  int m = C.size1();
  int n = C.size2();
  int k = A.size2();
  BOOST_UBLAS_CHECK (C.size1() == A.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (A.size2() == B.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (B.size2() == C.size2(), ublas::bad_size ());
  cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans,
	      m, n, k, 1.0, &A.data()[0], A.size2(), &B.data()[0], B.size2(), 1.0, &C.data()[0], C.size1());
}

BOOST_UBLAS_INLINE
void
cblas_prod (ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& C,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& A,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& B) {
  // C(m,n) <- aplha A(m,k) B(k,n) + beta C(m,n)
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  int m = C.size1();
  int n = C.size2();
  int k = A.size2();
  BOOST_UBLAS_CHECK (C.size1() == A.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (A.size2() == B.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (B.size2() == C.size2(), ublas::bad_size ());
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
	      m, n, k, 1.0, &A.data()[0], A.size1(), &B.data()[0], B.size2(), 1.0, &C.data()[0], C.size1());
}

BOOST_UBLAS_INLINE
void
cblas_prod (ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& C,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& A,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& B) {
  // C(m,n) <- aplha A(m,k) B(k,n) + beta C(m,n)
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  int m = C.size1();
  int n = C.size2();
  int k = A.size2();
  BOOST_UBLAS_CHECK (C.size1() == A.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (A.size2() == B.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (B.size2() == C.size2(), ublas::bad_size ());
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	      m, n, k, 1.0, &A.data()[0], A.size2(), &B.data()[0], B.size1(), 1.0, &C.data()[0], C.size1());
}

BOOST_UBLAS_INLINE
void
cblas_prod (ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& C,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& A,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& B) {
  // C(m,n) <- aplha A(m,k) B(k,n) + beta C(m,n)
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  int m = C.size1();
  int n = C.size2();
  int k = A.size2();
  BOOST_UBLAS_CHECK (C.size1() == A.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (A.size2() == B.size1(), ublas::bad_size ());
  BOOST_UBLAS_CHECK (B.size2() == C.size2(), ublas::bad_size ());
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
	      m, n, k, 1.0, &A.data()[0], A.size1(), &B.data()[0], B.size1(), 1.0, &C.data()[0], C.size1());
}



BOOST_UBLAS_INLINE
void
cblas_prodv (ublas::vector<double, ublas::unbounded_array<double> >& c,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& A,
	    const ublas::vector<double, ublas::unbounded_array<double> >& b) {
  //returns c = Ab + c
  int m = A.size1();
  int n = A.size2();
  cblas_dgemv(CblasColMajor, CblasNoTrans,
	      m, n, 1.0, &A.data()[0], A.size1(), &b.data()[0], 1, 1.0, &c.data()[0], 1);

}

BOOST_UBLAS_INLINE
void
cblas_prodv (ublas::vector<double, ublas::unbounded_array<double> >& c,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& A,
	    const ublas::vector<double, ublas::unbounded_array<double> >& b) {
  //returns c = Ab + c
  int m = A.size1();
  int n = A.size2();
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
	      m, n,
	      1.0,
	      &A.data()[0], A.size2(),
	      &b.data()[0], 1,
	      1.0,
	      &c.data()[0], 1);

}

BOOST_UBLAS_INLINE
void
cblas_prodv_trans (ublas::vector<double, ublas::unbounded_array<double> >& c,
	    const ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& A,
	    const ublas::vector<double, ublas::unbounded_array<double> >& b) {
  //returns c = Ab + c
  int m = A.size1();
  int n = A.size2();
  cblas_dgemv(CblasColMajor, CblasTrans,
	      m, n, 1.0, &A.data()[0], A.size1(), &b.data()[0], 1, 1.0, &c.data()[0], 1);

}

BOOST_UBLAS_INLINE
void
cblas_prodv_trans (ublas::vector<double, ublas::unbounded_array<double> >& c,
	    const ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& A,
	    const ublas::vector<double, ublas::unbounded_array<double> >& b) {
  //returns c = Ab + c
  int m = A.size1();
  int n = A.size2();
  cblas_dgemv(CblasRowMajor, CblasTrans,
	      m, n,
	      1.0,
	      &A.data()[0], A.size2(),
	      &b.data()[0], 1,
	      1.0,
	      &c.data()[0], 1);

}

BOOST_UBLAS_INLINE
void
cblas_outer_prod (ublas::matrix<double, ublas::column_major, ublas::unbounded_array<double> >& A,
		   const ublas::vector<double, ublas::unbounded_array<double> >& x,
		   const ublas::vector<double, ublas::unbounded_array<double> >& y){
  //returns A = x*y' + A;
  cblas_dger(CblasColMajor, A.size1(), A.size2(), 1, &x.data()[0], 1, &y.data()[0], 1, &A.data()[0], A.size1());
}

BOOST_UBLAS_INLINE
void
cblas_outer_prod (ublas::matrix<double, ublas::row_major, ublas::unbounded_array<double> >& A,
		   const ublas::vector<double, ublas::unbounded_array<double> >& x,
		   const ublas::vector<double, ublas::unbounded_array<double> >& y){
  //returns A = x*y' + A;
  cblas_dger(CblasRowMajor, A.size1(), A.size2(), 1, &x.data()[0], 1, &y.data()[0], 1, &A.data()[0], A.size2());
}

BOOST_UBLAS_INLINE
double
cblas_inner_prod(ublas::vector<double, ublas::unbounded_array<double> >& u1,
			ublas::vector<double, ublas::unbounded_array<double> >& u2){
  return cblas_ddot(u1.size(), &u1.data()[0], 1, &u2.data()[0], 1);
}




#endif /* CBLAS_MM_HPP_ */
