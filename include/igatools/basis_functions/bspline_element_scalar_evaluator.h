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
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/utils/dynamic_multi_array.h>


IGA_NAMESPACE_OPEN



/**
 * @brief blabla
 * @todo Document this class
 */
class Values1DConstView
{
public:
    /** Type for the container of one dimensional values on a single interval for a single scalar function.*/
    using Values1D = typename DenseMatrix::MatrixRowType ;

    using const_iterator = typename Values1D::const_iterator;

    Values1DConstView() = default;

    Values1DConstView(const DenseMatrix &funcs,const Index func_id)
        :
        funcs_(&funcs),
        func_id_(func_id)
    {
        Assert(func_id >= 0 && func_id < Size(funcs_->size1()),
               ExcIndexRange(func_id,0,Size(funcs_->size1())))
    }


    Values1DConstView(const Values1DConstView &view) = default ;
    Values1DConstView(Values1DConstView &&view) = default ;

    Values1DConstView &operator=(const Values1DConstView &view) = default;
    Values1DConstView &operator=(Values1DConstView &&view) = default;

    Real operator()(const Index point_id) const;

    Size get_num_points() const;

private:
    const DenseMatrix *funcs_ = nullptr;
    Index func_id_;
};



/**
 * @brief blabla
 * @todo Document this class
 */
template <int dim>
class
    BSplineElementScalarEvaluator
{
public:
    /** Type for the one dimensional values on a single interval for a single scalar function.*/
    using Values1D = typename DenseMatrix::MatrixRowType ;

    /**
     * Typedef for specifying the derivatives of the scalar basis function in the
     * reference domain.
     */
    template <int deriv_order>
    using Derivative = Derivatives<dim,1,1,deriv_order>;


    /** @name Constructors */
    ///@{
    BSplineElementScalarEvaluator() = delete;

    BSplineElementScalarEvaluator(const std::vector<std::array<Values1DConstView,dim>> &values1D);

    BSplineElementScalarEvaluator(const BSplineElementScalarEvaluator<dim> &bspline) = default;
    BSplineElementScalarEvaluator(BSplineElementScalarEvaluator<dim> &&bspline) = default;

    ~BSplineElementScalarEvaluator() = default;
    ///@}



    /** @name Assignment operators */
    ///@{
    BSplineElementScalarEvaluator<dim> &operator=(const BSplineElementScalarEvaluator<dim> &bspline) = default;
    BSplineElementScalarEvaluator<dim> &operator=(BSplineElementScalarEvaluator<dim> &&bspline) = default;
    ///@}


    Real evaluate_derivative(
        const TensorIndex<dim> &order_tensor_id,
        const TensorIndex<dim> &point_tensor_id) const;


    void evaluate_derivative_at_points(
        const TensorIndex<dim> &order_tensor_id,
        DynamicMultiArray<Real,dim> & derivatives) const;


    /** Returns the number of points in each direction for which the 1D values are associated. */
    TensorSize<dim> get_num_points() const;

    const std::array<Values1DConstView,dim>& get_derivative_components_view(const int order) const;

    const Values1DConstView & get_values_view(const int order,const int dir) const;

private:

    /**
     * values[i][j] are the values at the n_qp evaluation points of the i-th derivative
     * along the j-th direction.
     */
    std::vector<std::array<Values1DConstView,dim>> values1D_;
};



IGA_NAMESPACE_CLOSE


#endif // BSPLINE_ELEMENT_SCALAR_EVALUATOR_H_
