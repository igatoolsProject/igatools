//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

#include <igatools/base/quadrature_tensor_product.h>
#include <igatools/base/exceptions.h>
#include <igatools/utils/multi_array_utils.h>

IGA_NAMESPACE_OPEN

template<int dim_>
QuadratureTensorProduct<dim_>::
QuadratureTensorProduct()
    :
    EvaluationPoints<dim_>()
{
    this->weights_have_tensor_product_struct_ = true;
}



template<int dim_>
QuadratureTensorProduct<dim_>::
QuadratureTensorProduct(
    const TensorSize<dim_> num_points,
    void (*compute_coords_and_weight_1d_in)
    (const int n_pts_id, vector<Real> &coords,vector<Real> &weights),
    const Real eps_scaling)
    :
    EvaluationPoints<dim_>(),
    compute_coords_and_weight_1d(compute_coords_and_weight_1d_in)
{
    this->weights_have_tensor_product_struct_ = true;

    Assert(compute_coords_and_weight_1d != nullptr,ExcNullPtr());

    Assert(eps_scaling >= Real(0.0) && eps_scaling < Real(0.5),
           ExcMessage("The scaling factor must be >= 0.0 and < 0.5"));

    DirectionArray coords;
    DirectionArray weights_1d;
    for (int i = 0; i < dim_; ++i)
    {
        const auto n_pts = num_points[i];

        coords[i].resize(n_pts);
        weights_1d[i].resize(n_pts);
        compute_coords_and_weight_1d(n_pts, coords[i], weights_1d[i]);

        if (eps_scaling > 0)
            for (int ip = 0; ip < n_pts; ++ip)
                coords[i][ip] = 0.5 +
                                (coords[i][ip] / 0.5 - 1.0) * (0.5 - eps_scaling);
    }


    const int n_pts_total = num_points.flat_size();
    ValueVector<Point> points(n_pts_total);

    const auto n_pts_w = MultiArrayUtils<dim_>::compute_weight(num_points);
    for (int pt_flat_id = 0 ; pt_flat_id < n_pts_total ; ++pt_flat_id)
    {
        const auto pt_tensor_id =
            MultiArrayUtils<dim_>::flat_to_tensor_index(pt_flat_id,n_pts_w);

        for (int i = 0 ; i < dim_ ; ++i)
            points[pt_flat_id][i] = coords[i][pt_tensor_id[i]];
    }
    this->reset_points_coordinates_and_weights(points,weights_1d);
}



template<int dim_>
QuadratureTensorProduct<dim_>::
QuadratureTensorProduct(
    const PointVector &points,
    const DirectionArray &weights_1d,
    const BBox<dim_> &bounding_box)
    :
    EvaluationPoints<dim_>(points,weights_1d,bounding_box)
{}

IGA_NAMESPACE_CLOSE

#include <igatools/base/quadrature_tensor_product.inst>
