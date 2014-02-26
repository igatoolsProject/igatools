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

using std::vector;
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
Quadrature(const TensorProductArray<dim> &points,
           const TensorProductArray<dim> &weights)
    :
    points_(points),
    weights_(weights)
{
    Assert(points.tensor_size() == weights.tensor_size(),
           ExcMessage("Wrong sizes."));
}




template< int dim >
TensorProductArray<dim>
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

    out << "weights (flat tensor product):" << endl;
    out << get_weights().get_flat_tensor_product() << endl;

    out << "coordinates:" << endl;
    this->get_points().print_info(out);
    out << endl;

    out << "coordinates (flat cartesian_product):" << endl;
    out << get_points().get_flat_cartesian_product() << endl;

    out << endl;
}


template< int dim >
Quadrature<dim>
Quadrature<dim> ::
get_restriction(const int face_id) const
{
    AssertIndexRange(face_id, UnitElement<dim>::faces_per_element);

    const int face_const_dir = UnitElement<dim>::face_constant_direction[face_id];
    const Real face_value = UnitElement<dim>::face_constant_coordinate[face_id];

    const auto points_old  = this->get_points();
    const auto weights_old = this->get_weights();

    TensorProductArray<dim> points_new;
    TensorProductArray<dim> weights_new;

    points_new.copy_data_direction(face_const_dir,vector<Real>(1,face_value));
    weights_new.copy_data_direction(face_const_dir,vector<Real>(1,1.0));

    for (auto i : UnitElement<dim>::face_active_directions[face_id])
    {
        points_new.copy_data_direction(i,points_old.get_data_direction(i));
        weights_new.copy_data_direction(i,weights_old.get_data_direction(i));
    }

    return Quadrature<dim>(points_new, weights_new);
}

template< int dim >
Quadrature< dim+1 > extend_quad_dim(const Quadrature< dim > &quad_surf,
                                    const int face_id)
{
    AssertIndexRange(face_id, UnitElement<dim+1>::faces_per_element);


    const int const_direction = UnitElement<dim+1>::face_constant_direction[face_id];

    TensorProductArray<dim> points_old = quad_surf.get_points();
    TensorProductArray<dim> weights_old = quad_surf.get_weights();

    TensorProductArray<dim+1> points_new;
    TensorProductArray<dim+1> weights_new;

    for (int i = 0, j = 0; i < dim + 1; ++i, ++j)
    {
        if (i != const_direction)
        {
            points_new.copy_data_direction(i,points_old.get_data_direction(j));
            weights_new.copy_data_direction(i,weights_old.get_data_direction(j));
        }
        else
        {
            points_new.copy_data_direction(
                i,vector<Real>(1,UnitElement<dim+1>::face_constant_coordinate[face_id]));
            weights_new.copy_data_direction(i,vector<Real>(1,1.0));
            --j;
        }
    }

    return Quadrature<dim+1>(points_new, weights_new);
}


#if 0
template< int dim >
Quadrature< dim > restricted_quad(const Quadrature< dim > &quad_surf,
                                  const int face_id)
{
    AssertIndexRange(face_id, UnitElement<dim>::faces_per_element);

    const int face_const_dir = UnitElement<dim>::face_constant_direction[face_id];
    const Real face_value = UnitElement<dim>::face_constant_coordinate[face_id];

    const auto points_old  = quad_surf.get_points();
    const auto weights_old = quad_surf.get_weights();

    TensorProductArray<dim> points_new;
    TensorProductArray<dim> weights_new;

    points_new.copy_data_direction(face_const_dir,vector<Real>(1,face_value));
    weights_new.copy_data_direction(face_const_dir,vector<Real>(1,1.0));

    for (auto i : UnitElement<dim>::face_active_directions[face_id])
    {
        points_new.copy_data_direction(i,points_old.get_data_direction(i));
        weights_new.copy_data_direction(i,weights_old.get_data_direction(i));
    }

    return Quadrature<dim>(points_new, weights_new);
}
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/base/quadrature.inst>
