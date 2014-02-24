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

#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>

#include <vector>

IGA_NAMESPACE_OPEN

/**
 *  This class represents a function
 *  \f[ f: \mathbb{R}^n \to \mathbb{R}^{\underbrace{m \times \dots \times m}_{r times}} \f].
 *
 *  For example:
 *  - Function<n,m,1> is an m-vector valued function from <i> R <sup> n </sup> </i>.
 *  - Function<n,m,0> is a scalar function for any m.
 *
 *  This is a pure virtual function and to be instantiate a derived class
 *  defining evaluate must be provided.
 *
 *  For example to define a liner function from R^m to R^n
 *  \code
 *  template<int m, int n>
 *  class LinearFunction : public Function<m, n, 1 >
 *  {
 *  public:
 *    LinearFunction();
 *    void evaluate(std::vector < typename LinearFunction<dim, rdim>::Point >  & points,
 *                  std::vector     < typename LinearFunction<dim, rdim>::ValueType >  & values
 *                 ) const
 *    {
 *      for (int i=0; i<points.size(); i++)
 *          values[i] = action(A,points[i]) + b;
 *    }
 *  private:
 *    Tensor <dim, 1, tensor::covariant, Tensor<rdim, 1, tensor::contravariant, Real> > A;
 *    Tensor<rdim, 1, tensor::contravariant, Real> b;
 *  };
 *  \endcode
 */

template<int dim, int dim_range = 1, int rank = 1>
class Function
{
public:
    /** Types for the input/output evaluation arguments */
    ///@{
    /**
     * Type for the input argument of the function.
     */
    using PointType = Point<dim>;

    /**
     * Type for the return of the function.
     */
    using ValueType = Values<dim, dim_range, rank>;

    /**
     * Type for the gradient of the function.
     */
    using GradientType = Derivatives<dim, dim_range, rank, 1>;

    /**
     * Type for the hessian of the function.
     */
    using HessianType = Derivatives<dim, dim_range, rank, 2>;
    ///@}

    /** @name Constructors and destructor. */
    ///@{
    /** Constructor */
    Function() = default;

    /** Destructor */
    virtual ~Function();
    ///@}

    /** @name Evaluations of the Function at some points */
    ///@{
    /** Compute the @p values of Function at some @p points. */
    virtual void evaluate(const std::vector<PointType> &points,
                          std::vector<ValueType> &values) const = 0;

    /** Compute the @p gradients of Function at some @p points. */
    virtual void evaluate_gradients(
        const std::vector<PointType> &points,
        std::vector<GradientType> &gradient) const;

    /** Compute the @p hessians of Function at some @p points. */
    virtual void evaluate_hessians(
        const std::vector<PointType> &points,
        std::vector<HessianType> &hessians) const;

    /** Compute the @values and the @p gradients of Function at some @p points. */
    virtual void evaluate_values_and_gradients(
        const std::vector<PointType> &points,
        std::vector<ValueType> &values,
        std::vector<GradientType> &gradients) const;
    ///@}
};


/**
 * Scalar Function.
 */
template<int dim>
using ScalarFunction = Function<dim, 1, 1>;

/**
 * Vector Function.
 */
template<int dim>
using VectorFunction = Function<dim, dim, 1>;

/**
 * Tensor Function.
 */
template<int dim>
using TensorFunction = Function<dim, dim, 2>;


IGA_NAMESPACE_CLOSE

#endif /* _FUNCTIONS_H */
