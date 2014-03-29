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

#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>
#include <igatools/base/exceptions.h>


IGA_NAMESPACE_OPEN

Size
Values1DConstView::
get_num_points() const
{
	return funcs_->size2();
}

Real
Values1DConstView::
operator()(const Index point_id) const
{
    Assert(point_id >= 0 && point_id < this->get_num_points(),
           ExcIndexRange(point_id,0,this->get_num_points()));

    return (*funcs_)(func_id_,point_id);
}


template <int dim>
BSplineElementScalarEvaluator<dim>::
BSplineElementScalarEvaluator(const std::vector<std::array<Values1DConstView,dim>> &values1D)
    :
    values1D_(values1D)
{
    Assert(!values1D_.empty(),ExcEmptyObject());
}


template <int dim>
const std::array<Values1DConstView,dim>&
BSplineElementScalarEvaluator<dim>::
get_derivative_components_view(const int order) const
{
	Assert(order >= 0 && order< values1D_.size(),
			ExcIndexRange(order,0,values1D_.size()));
	return values1D_[order];
}

template <int dim>
const Values1DConstView &
BSplineElementScalarEvaluator<dim>::
get_values_view(const int order,const int dir) const
{
	Assert(dir >= 0 && dir < dim,
			ExcIndexRange(dir,0,dim));
	return get_derivative_components_view(order)[dim];
}




template <int dim>
TensorSize<dim>
BSplineElementScalarEvaluator<dim>::
get_num_points() const
{
	Assert(!values1D_.empty(),ExcEmptyObject());

	TensorSize<dim> n_points;
	for (int i = 0 ; i < dim ; ++i)
	{
		n_points(i) = values1D_[0][i].get_num_points();

#ifndef NDEBUG
		for (const auto &values : values1D_)
			Assert(n_points(i) == values[i].get_num_points(),
					ExcDimensionMismatch(n_points(i),values[i].get_num_points()));
#endif
	}
	return n_points;
}


template <int dim>
Real
BSplineElementScalarEvaluator<dim>::
evaluate_derivative(
    const TensorIndex<dim> &order_tensor_id,
    const TensorIndex<dim> &point_tensor_id) const
{
#ifndef NDEBUG
    for (int i = 0 ; i < dim ; ++i)
        Assert(order_tensor_id[i] >= 0 && order_tensor_id[i] < values1D_.size(),
               ExcIndexRange(order_tensor_id[i],0,values1D_.size()));
#endif

    Assert(dim > 0, ExcLowerRange(dim,1));

    // Main computation
    Real partial_derivative = values1D_[order_tensor_id[0]][0](point_tensor_id[0]);
    for (int i = 1; i < dim; ++i)
    {
        partial_derivative *= values1D_[order_tensor_id[i]][i](point_tensor_id[i]);
    }


    return partial_derivative;
}

template <>
Real
BSplineElementScalarEvaluator<0>::
evaluate_derivative(
    const TensorIndex<0> &order_tensor_id,
    const TensorIndex<0> &point_tensor_id) const
{
    return 0.0;
}

template <int dim>
void
BSplineElementScalarEvaluator<dim>::
evaluate_derivative_at_points(
    const TensorIndex<dim> &order_tensor_id,
    DynamicMultiArray<Real,dim> & derivatives) const
{
	TensorSize<dim> n_points = this->get_num_points();


	Assert(false,ExcNotImplemented());
	AssertThrow(false,ExcNotImplemented());
}


IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/bspline_element_scalar_evaluator.inst>
