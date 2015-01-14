/*
 * boostHeader.h
 *
 *  Created on: Mar 13, 2010
 *      Author: bhole
 */

#ifndef BOOSTHEADER_H_
#define BOOSTHEADER_H_


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>


#define vectorB boost::numeric::ublas::vector
#define matrixB boost::numeric::ublas::matrix

#define matrix_rowB boost::numeric::ublas::matrix_row
#define matrix_columnB boost::numeric::ublas::matrix_column

#define column_majorB boost::numeric::ublas::column_major


#endif /* BOOSTHEADER_H_ */
