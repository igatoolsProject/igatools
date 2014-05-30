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


template <int dim>
BSplineElementScalarEvaluator<dim>::
BSplineElementScalarEvaluator(const std::vector<std::array<Values1DConstView,dim>> &values1D)
    :
    values1D_(values1D)
{
    Assert(!values1D_.empty(),ExcEmptyObject());
}


template <int dim>
const std::array<Values1DConstView,dim> &
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
/*
template <int dim>
template <int k>
void
BSplineElementScalarEvaluator<dim>::
recursive_multiplication(
    const TensorIndex<dim> &order_tensor_id,
    DynamicMultiArray<Real,dim-k+1> & derivative) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());

    const int old_rank = dim-k+1;
    const int new_dir = k-1;

    TensorSize<rank> num_pts;
    TensorSize<rank+1> num_pts_new;
    for ( int rank = 0 ; rank < old_rank ; ++rank)
    {
        dir = dim-rank-1;
        num_pts(rank) = this->get_num_points()(dir);
        num_pts_new(rank+1) = num_pts(rank);

        Assert(num_pts(rank) == derivative.tensor_size()(dir),
                ExcDimensionMismatch(num_pts(rank),derivative.tensor_size()(dir)));
    }
    const Size pts_flat_size = num_pts.flat_size();

    TensorIndex<rank> pt_tensor_id;
    TensorIndex<rank> pt_tensor_w = MultiArrayUtils<rank>::compute_weight(num_pts);

    //number of points in the new direction
    const Size n_pts_new_dir = this->get_num_points()(new_dir);
    num_pts_new(0) = n_pts_new_dir;

    const auto &deriv_new_dir = values1D_[order_tensor_id[new_dir]][new_dir];
    Index pt_flat_id_new = 0;
    for (Index pt_flat_id = 0 ; pt_flat_id < pts_flat_size ; ++pt_flat_id)
    {
        //loop over the previous points
        const Real value = derivative(pt_flat_id);

        for (Index i = 0 ; i < n_pts_new_dir ; ++i,++pt_flat_id_new)
        {
            //loop over the new point index
            dervative_new(pt_flat_id_new) = value * deriv_new_dir(i);
        }

    }

    recursive_multiplication<k+1>(order_tensor_id,derivative_next);
}
//*/

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


template <>
void
BSplineElementScalarEvaluator<1>::
evaluate_derivative_at_points(
    const TensorIndex<1> &order_tensor_id,
    DynamicMultiArray<Real,1> &derivatives) const
{
    TensorSize<1> n_points = this->get_num_points();

    Assert(derivatives.tensor_size()(0) == n_points(0),
           ExcDimensionMismatch(derivatives.tensor_size()(0),n_points(0)));

    const auto &deriv_dir_0 = this->get_derivative_components_view(order_tensor_id[0])[0];
    for (Index flat_pt_id_0 = 0 ; flat_pt_id_0 < n_points(0) ; ++flat_pt_id_0)
        derivatives(flat_pt_id_0) = deriv_dir_0(flat_pt_id_0);
}

template <>
void
BSplineElementScalarEvaluator<2>::
evaluate_derivative_at_points(
    const TensorIndex<2> &order_tensor_id,
    DynamicMultiArray<Real,2> &derivatives) const
{
    TensorSize<2> n_points = this->get_num_points();

#ifndef NDEBUG
    for (int dir = 0 ; dir < 2 ; ++dir)
        Assert(derivatives.tensor_size()(dir) == n_points(dir),
               ExcDimensionMismatch(derivatives.tensor_size()(dir),n_points(dir)));
#endif

    const auto &deriv_dir_1 = this->get_derivative_components_view(order_tensor_id[1])[1];
    const auto &deriv_dir_0 = this->get_derivative_components_view(order_tensor_id[0])[0];

    Index flat_pt_id = 0;
    for (Index flat_pt_id_1 = 0 ; flat_pt_id_1 < n_points(1) ; ++flat_pt_id_1)
    {
        const Real deriv_dir_1_pt = deriv_dir_1(flat_pt_id_1);

        for (Index flat_pt_id_0 = 0 ; flat_pt_id_0 < n_points(0) ; ++flat_pt_id_0)
            derivatives(flat_pt_id++) = deriv_dir_0(flat_pt_id_0) * deriv_dir_1_pt;
    }
}

template <>
void
BSplineElementScalarEvaluator<3>::
evaluate_derivative_at_points(
    const TensorIndex<3> &order_tensor_id,
    DynamicMultiArray<Real,3> &derivatives) const
{
    TensorSize<3> n_points = this->get_num_points();

#ifndef NDEBUG
    for (int dir = 0 ; dir < 3 ; ++dir)
        Assert(derivatives.tensor_size()(dir) == n_points(dir),
               ExcDimensionMismatch(derivatives.tensor_size()(dir),n_points(dir)));
#endif

    const auto &deriv_dir_2 = this->get_derivative_components_view(order_tensor_id[2])[2];
    const auto &deriv_dir_1 = this->get_derivative_components_view(order_tensor_id[1])[1];
    const auto &deriv_dir_0 = this->get_derivative_components_view(order_tensor_id[0])[0];

    Index flat_pt_id = 0;
    for (Index flat_pt_id_2 = 0 ; flat_pt_id_2 < n_points(2) ; ++flat_pt_id_2)
    {
        const Real deriv_dir_2_pt = deriv_dir_2(flat_pt_id_2);

        for (Index flat_pt_id_1 = 0 ; flat_pt_id_1 < n_points(1) ; ++flat_pt_id_1)
        {
            const Real deriv_dir_1_pt =  deriv_dir_1(flat_pt_id_1);
            const Real deriv_old_dirs_pt = deriv_dir_2_pt * deriv_dir_1_pt;

            for (Index flat_pt_id_0 = 0 ; flat_pt_id_0 < n_points(0) ; ++flat_pt_id_0)
                derivatives(flat_pt_id++) = deriv_dir_0(flat_pt_id_0) * deriv_old_dirs_pt;
        }
    }
}


template <int dim>
void
BSplineElementScalarEvaluator<dim>::
evaluate_derivative_at_points(
    const TensorIndex<dim> &order_tensor_id,
    DynamicMultiArray<Real,dim> &derivatives) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}


IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/bspline_element_scalar_evaluator.inst>
