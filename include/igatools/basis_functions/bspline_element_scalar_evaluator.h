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


#ifndef BSPLINE_ELEMENT_SCALAR_EVALUATOR_H_
#define BSPLINE_ELEMENT_SCALAR_EVALUATOR_H_


#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/basis_functions/values1d_const_view.h>
#include <igatools/utils/dynamic_multi_array.h>


IGA_NAMESPACE_OPEN


/**
 * @brief Evaluator of scalar BSpline values and derivatives from
 * one-dimensional values.
 *
 * This purpose of this class is to isolate the methods for computing
 * the BSpline values and derivatives
 * from the (precomputed) one-dimensional values (and derivative).
 * Each instance of this class represents a single
 * BSpline \f$ N_i : \mathbb{R}^{dim} \to \mathbb{R}\f$ that
 * must be evaluated at some points arranged in a tensor-product way.
 *
 * @todo Document more
 * @author M. Martinelli
 * @date 29 Mar 2014
 */
template <int dim>
class BSplineElementScalarEvaluator
{
public:
    /** Type for the one dimensional values on a single interval for
     * a single scalar function.*/
    using Values1D = typename DenseMatrix::MatrixRowType ;

    /**
     * Typedef for specifying the derivatives of the scalar basis function in the
     * reference domain.
     */
    template <int order>
    using Derivative = Derivatives<dim,1,1,order>;


    /** @name Constructors */
    ///@{
    /** Default constructor. Not allowed to be used. */
    BSplineElementScalarEvaluator() = delete;

    /**
     * Constructor.
     * It builds the scalar BSpline evaluator from the one-dimensional views of values and derivatives.
     * @p values1D[i][j] are the <tt>i</tt>-th order (one-dimensional) derivatives along the <tt>j</tt>-th direction.
     */
    BSplineElementScalarEvaluator(const vector<std::array<Values1DConstView,dim>> &values1D);

    /** Copy constructor. */
    BSplineElementScalarEvaluator(const BSplineElementScalarEvaluator<dim> &bspline) = default;

    /** Move constructor. */
    BSplineElementScalarEvaluator(BSplineElementScalarEvaluator<dim> &&bspline) = default;

    /** Destructor. */
    ~BSplineElementScalarEvaluator() = default;
    ///@}



    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator. */
    BSplineElementScalarEvaluator<dim> &operator=(const BSplineElementScalarEvaluator<dim> &bspline) = default;

    /** Move assignment operator. */
    BSplineElementScalarEvaluator<dim> &operator=(BSplineElementScalarEvaluator<dim> &&bspline) = default;
    ///@}


    /**
     * Evaluate and returns one partial derivative in one point.
     * The order of the partial derivative is specified by the tensor-index @p order_tensor_id,
     * while the point is specified by its tensor-index @p point_tensor_id.
     */
    Real evaluate_derivative(
        const TensorIndex<dim> &order_tensor_id,
        const TensorIndex<dim> &point_tensor_id) const;


    /**
     * Evaluate and returns one partial @p derivative in all points.
     * The order of the partial derivative is specified by the tensor-index @p order_tensor_id.
     */
    void evaluate_derivative_at_points(
        const TensorIndex<dim> &order_tensor_id,
        DynamicMultiArray<Real,dim> &derivative) const;


    /** Returns the number of points in each direction for which the 1D values are associated. */
    TensorSize<dim> get_num_points() const;

    /**
     * Returns a const view to the one-dimensional derivatives of given @p order,
     * along all coordinate directions.
     */
    const std::array<Values1DConstView,dim> &get_derivative_components_view(const int order) const;

    /**
     * Returns a const view to the one-dimensional derivative of given @p order,
     * along the coordinate direction @p dir.
     */
    const Values1DConstView &get_values_view(const int order,const int dir) const;

private:
    /*
        template <int k>
        void recursive_multiplication(
            const TensorIndex<dim> &order_tensor_id,
            DynamicMultiArray<Real,dim> & derivative) const;
    //*/
    /**
     * values[i][j] are the values at the n_qp evaluation points of the i-th derivative
     * along the j-th direction.
     */
    vector<std::array<Values1DConstView,dim>> values1D_;
};



IGA_NAMESPACE_CLOSE


#endif // BSPLINE_ELEMENT_SCALAR_EVALUATOR_H_
