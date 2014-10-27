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

#include <igatools/base/quadrature.h>
#include <igatools/base/exceptions.h>
#include <igatools/geometry/unit_element.h>

using std::array;
using std::endl;

IGA_NAMESPACE_OPEN

template< int dim >
Quadrature< dim >::Quadrature(const TensorSize<dim> num_points)
    :
    points_(num_points),
    weights_(num_points)
{}



template< int dim >
Quadrature< dim >::Quadrature(const Index num_points)
    :
    points_(num_points),
    weights_(num_points)
{}



template< int dim >
Quadrature< dim >::
Quadrature(const CartesianProductArray<Real,dim> &points,
           const TensorProductArray<dim> &weights)
    :
    points_(points),
    weights_(weights)
{
    Assert(points.tensor_size() == weights.tensor_size(),
           ExcMessage("Sizes of points and weights do not match."));
}




template< int dim >
CartesianProductArray<Real,dim>
Quadrature< dim >::
get_points() const noexcept
{
    return points_;
}



template< int dim >
TensorProductArray<dim>
Quadrature< dim >::
get_weights() const noexcept
{
    return weights_;
}



template< int dim >
TensorSize<dim>
Quadrature< dim >::
get_num_points_direction() const noexcept
{
    return points_.tensor_size();
}



template< int dim >
Size
Quadrature< dim >::
get_num_points() const noexcept
{
    return points_.flat_size();
}



template< int dim >
void
Quadrature< dim >::
print_info(LogStream &out) const
{
    out << "Number of points:" << get_num_points() << endl;

    out << "weights:" << endl;
    get_weights().print_info(out);
    out << endl;

    // TODO (pauletti, Aug 26, 2014): redundant info, remove
    out << "weights (flat tensor product):" << endl;
    get_weights().get_flat_tensor_product().print_info(out);
    out << endl;

    out << "coordinates:" << endl;
    this->get_points().print_info(out);
    out << endl;

    // TODO (pauletti, Aug 26, 2014): redundant info, remove
    out << "coordinates (flat cartesian_product):" << endl;
    get_points().get_flat_cartesian_product().print_info(out);
    out << endl;

    out << endl;
}

template< int dim >
template<int k>
Quadrature<dim>
Quadrature<dim>::
collapse_to_sub_element(const int sub_elem_id) const
{
    auto &k_elem = UnitElement<dim>::template get_elem<k>(sub_elem_id);

    PointArray  new_points;
    WeigthArray new_weights;

    const int n_dir = k_elem.constant_directions.size();
    for (int j=0; j<n_dir; ++j)
    {
        auto dir = k_elem.constant_directions[j];
        auto val = k_elem.constant_values[j];
        new_points.copy_data_direction(dir,vector<Real>(1, val));
        new_weights.copy_data_direction(dir,vector<Real>(1, 1.0));

    }

    for (auto i : k_elem.active_directions)
    {
        new_points.copy_data_direction(i,points_.get_data_direction(i));
        new_weights.copy_data_direction(i,weights_.get_data_direction(i));
    }

    return Quadrature<dim>(new_points, new_weights);
}



#if 0
template< int face_dim >
Quadrature<face_dim+1> extend_face_quad(const Quadrature<face_dim> &quad,
                                        const int face_id)
{
    const int dim = face_dim + 1;

    AssertIndexRange(face_id, UnitElement<dim>::faces_per_element);
    const int const_direction = UnitElement<dim>::face_constant_direction[face_id];

    CartesianProductArray<Real,face_dim> points_old = quad.get_points();
    TensorProductArray<face_dim> weights_old = quad.get_weights();

    CartesianProductArray<Real,dim> points_new;
    TensorProductArray<dim> weights_new;

    for (int i = 0, j = 0; i < dim ; ++i, ++j)
    {
        if (i != const_direction)
        {
            points_new.copy_data_direction(i,points_old.get_data_direction(j));
            weights_new.copy_data_direction(i,weights_old.get_data_direction(j));
        }
        else
        {
            points_new.copy_data_direction(
                i,
                vector<Real>(1,UnitElement<dim>::face_constant_coordinate[face_id]));
            weights_new.copy_data_direction(i,vector<Real>(1,1.0));
            --j;
        }
    }

    return Quadrature<dim>(points_new, weights_new);
}
#endif
IGA_NAMESPACE_CLOSE

#include <igatools/base/quadrature.inst>
