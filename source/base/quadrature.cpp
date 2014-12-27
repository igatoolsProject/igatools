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

#include <igatools/base/quadrature.h>
#include <igatools/base/exceptions.h>
#include <igatools/geometry/unit_element.h>

using std::array;
using std::endl;

IGA_NAMESPACE_OPEN

template<int dim_>
Quadrature<dim_>::Quadrature(const TensorSize<dim> num_points)
    :
    points_(num_points),
    weights_(num_points)
{}



template<int dim_>
Quadrature<dim_>::Quadrature(const Index num_points)
    :
    points_(num_points),
    weights_(num_points)
{}



template<int dim_>
Quadrature<dim_>::
Quadrature(const CartesianProductArray<Real,dim> &points,
           const TensorProductArray<dim> &weights)
    :
    points_(points),
    weights_(weights)
{
    Assert(points.tensor_size() == weights.tensor_size(),
           ExcMessage("Sizes of points and weights do not match."));
}




template<int dim_>
auto
Quadrature<dim_>::get_points() const noexcept -> PointArray
{
    return points_;
}



template<int dim_>
auto
Quadrature<dim_>::get_weights() const noexcept -> WeigthArray
{
    return weights_;
}



template<int dim_>
auto
Quadrature<dim_>::
get_num_points_direction() const noexcept -> TensorSize<dim>
{
    return points_.tensor_size();
}



template<int dim_>
Size
Quadrature<dim_>::
get_num_points() const noexcept
{
    return points_.flat_size();
}



template<int dim_>
void
Quadrature<dim_>::
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

template<int dim_>
template<int k>
auto
Quadrature<dim_>::
collapse_to_sub_element(const int sub_elem_id) const -> self_t
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



template<int k, int dim>
Quadrature<dim>
extend_sub_elem_quad(const Quadrature<k> &quad,
                     const int sub_elem_id)
{

    auto &k_elem = UnitElement<dim>::template get_elem<k>(sub_elem_id);

    typename Quadrature<dim>::PointArray  new_points;
    typename Quadrature<dim>::WeigthArray new_weights;

    const int n_dir = k_elem.constant_directions.size();
    for (int j=0; j<n_dir; ++j)
    {
        auto dir = k_elem.constant_directions[j];
        auto val = k_elem.constant_values[j];
        new_points.copy_data_direction(dir,vector<Real>(1, val));
        new_weights.copy_data_direction(dir,vector<Real>(1, 1.0));

    }

    int ind = 0;
    for (auto i : k_elem.active_directions)
    {
        new_points.copy_data_direction(i,quad.get_points().get_data_direction(ind));
        new_weights.copy_data_direction(i,quad.get_weights().get_data_direction(ind));
        ++ind;
    }

    return Quadrature<dim>(new_points, new_weights);
}

IGA_NAMESPACE_CLOSE

#include <igatools/base/quadrature.inst>
