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
#include <igatools/utils/value_vector.h>

#include <vector>

IGA_NAMESPACE_OPEN

/**
 *  @brief Scalar, vector or tensor valued function
 *
 *  Function
 *  \f$ f: \mathbb{R}^n \to \mathbb{R}^{\underbrace{m \times \dots \times m}_{r_{times}}}. \f$
 *
 *  For example:
 *  - Function<n> is a scalar valued function on <i> R <sup> n </sup> </i>.
 *  - Function<n,m> is an m-dimensional vector valued function on
 *    <i> R <sup> n </sup> </i>.
 *  - Function<n,m,2> is a tensor valued function on <i> R <sup> n </sup> </i>.
 *
 *  This is a pure abstract class and cannot be instantiated.
 *  Its purpose is to give an unified interface for concrete classes implementing
 *  function evaluation through the specialization of Function::evaluate().
 *
 *  For example to define a linear function from R^m to R^n
 *  \code
 *  template<int m, int n>
 *  class LinearFunction : public Function<m, n>
 *  {
 *  public:
 *  using Function<m, n>::
 *    LinearFunction();
 *    void evaluate(ValueVector< typename LinearFunction<dim,rdim>::Point >  & points,
 *                  ValueVector< typename LinearFunction<dim,rdim>::Value >  & values
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

template<int dim, int range = 1, int rank = 1>
class Function
{
public:
    /** Types for the input/output evaluation arguments */
    ///@{
    /**
     * Type for the input argument of the function.
     */
    using Point = Points<dim>;

    /**
     * Type for the return of the function.
     */
    using Value = Values<dim, range, rank>;

    /**
     * Type for the derivative of the function.
     */
    template <int order>
    using Derivative = Derivatives<dim, range, rank, order>;

    /**
     * Type for the gradient of the function.
     */
    using Gradient = Derivative<1>;

    /**
     * Type for the hessian of the function.
     */
    using Hessian = Derivative<2>;

    /**
     * Typedef for specifying the divergence of the basis function.
     */
    using Div = Values<dim, range, rank-1>;
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
    virtual void evaluate(const ValueVector<Point> &points,
                          ValueVector<Value> &values) const = 0;

    /** Compute the @p gradients of Function at some @p points. */
    virtual void evaluate_gradients(const ValueVector<Point> &points,
                                    ValueVector<Gradient> &gradient) const;

    /** Compute the @p hessians of Function at some @p points. */
    virtual void evaluate_hessians(const ValueVector<Point> &points,
                                   ValueVector<Hessian> &hessians) const;

    /** Compute the @p values and the @p gradients of Function at some
     *  @p points. */
    virtual void evaluate_values_and_gradients(const ValueVector<Point> &points,
                                               ValueVector<Value> &values,
                                               ValueVector<Gradient> &gradients) const;
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

#endif /* __FUNCTIONS_H */
