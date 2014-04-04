//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
//
// This file is part of the igatools library.
//
// The igatools library is free software: you can use it, redistribute
// it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-+--------------------------------------------------------------------


#ifndef DENSE_MATRIX_H_
#define DENSE_MATRIX_H_

#include <igatools/base/config.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


IGA_NAMESPACE_OPEN

/** @addtogroup linear_algebra
 *@{
 */
/**
 * @brief Dense matrix with real entries.
 *
 * A typical use of it is the local matrix.
 * @note This class inherits (also the constructors) from the matrix type in the
 * <a href="http://www.boost.org/doc/libs/1_55_0/libs/numeric/ublas/doc/index.htm">Boost uBlas library</a>
 * Therefore it can be used as a boost::numeric::ublas::matrix<Real>.
 *
 * @author M. Martinelli, 2012, 2013, 2014
 * @author S. Pauletti, 2012, 2013
 */
class DenseMatrix : public boost::numeric::ublas::matrix<Real>
{
public:
    /** Type of the base class. */
    using BoostMatrix =  boost::numeric::ublas::matrix<Real> ;

    /** We inherith the constructors of the base class. */
    using BoostMatrix::BoostMatrix;

    /**
     * Assignment operator for assigning zeros to all entries of the matrix
     * by writing
     * @code
       DenseMatrix matrix;
       ... // working with the matrix
       matrix = 0.0; // reset the matrix entries to zero
       @endcode
     * @note If used in Debug mode with a @p value different from zero,
     * an assertion will be raised.
     */
    DenseMatrix &operator=(const Real value) ;

    /** Type for a row of the matrix. */
    using MatrixRowType = boost::numeric::ublas::matrix_row<const BoostMatrix> ;

    /** Returns the mtrix row identified by @p row_id. */
    MatrixRowType
    get_row(const int row_id) const;


    /**
     * Matrix inversion routine.
     * @note Uses lu_factorize and lu_substitute in uBLAS to invert a matrix
     */
    DenseMatrix inverse();


    /**
     * Returns the Frobenius norm of the matrix.
     */
    Real norm_frobenius() const;


    /** Returns the number of rows of the matrix. */
    Size get_num_rows() const;

    /** Returns the number of columns of the matrix. */
    Size get_num_cols() const;
};

/**@}*/

IGA_NAMESPACE_CLOSE

#endif /* DENSE_MATRIX_H_ */
