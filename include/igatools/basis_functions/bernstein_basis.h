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

#ifndef __BERNSTEIN_BASIS_H
#define __BERNSTEIN_BASIS_H

#include <igatools/base/config.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <igatools/utils/vector.h>

//TODO: document the implementation and expand introduction

IGA_NAMESPACE_OPEN


/**
 * The Bernstein polynomial basis of degree p.
 * The i-th (0<=i<=p) Bernstein polynomial
 * is defined as: \f$ B^p_i(x) = {p \choose i } x^i (1-x)^{n-i} \f$.
 *
 * They are defined in [0,1] and in the context of IGA
 * are used to generate B-Splines and NURBS through
 * the Bezier extraction.
 *
 * @note the return type for values and derivatives is
 * a matrix because of it expected use to generate the
 * B-spline. A typical example code of its use is:
 * \code
 * todo
 * \endcode
 *
 * @author pauletti 2013
 * @author M. Martinelli 2014
 *
 */
namespace BernsteinBasis
{
/**
 * Compute the values of all the p+1 Bernstein basis
 * of degree @p p at the given point \f$ x \f$.
 *
 * The values are returned in a vector B, where
 *  B[i] is \f$ B^p_i(x) \f$.
 *
 * @warning The point \f$ x \f$ must belong to the unit interval [0,1], otherwise an
 * assertion will be raised in Debug mode.
 */
boost::numeric::ublas::vector<Real>
evaluate(const int p, const Real x) ;

/**
 * Compute the values of all the p+1 Bernstein basis
 * of degree @p p at the given points \f$ x \f$.
 *
 * The values are returned in a matrix B, where
 *  B[i][j] is \f$ B^p_i(x_j) \f$.
 *
 * @warning The points \f$ x \f$ must belong to the unit interval [0,1], otherwise an
 * assertion will be raised in Debug mode.
 *
 */
boost::numeric::ublas::matrix<Real>
evaluate(const int p, const vector<Real> &x) ;


/**
 * Computes the derivatives of order k of the Berstein basis
 * at the given point \f$ x \f$.
 * The values are returned in a vector DkB, where
 *  DkB[i] is \f$ \frac{d^k B^p_i(x)}{dx^k} \f$.
 *
 * @warning The point \f$ x \f$ must belong to the unit interval [0,1], otherwise an
 * assertion will be raised in Debug mode.
 */
boost::numeric::ublas::vector<Real>
derivative(int k, const int p, const Real x) ;

/**
 * Computes the derivatives of order k of the Berstein basis
 * at the given points \f$ x \f$.
 * The values are returned in a matrix DkB, where
 *  DkB[i][j] is \f$ \frac{d^k B^p_i(x_j)}{dx^k} \f$.
 *
 * @warning The points \f$ x \f$ must belong to the unit interval [0,1], otherwise an
 * assertion will be raised in Debug mode.
 *
 */
boost::numeric::ublas::matrix<Real>
derivative(int k, const int p, const vector<Real> &x) ;


};


IGA_NAMESPACE_CLOSE


#endif
