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

#include <igatools/geometry/grid_uniform_quad_cache.h>


using std::shared_ptr;

IGA_NAMESPACE_OPEN

template <int dim_>
GridUniformQuadCache<dim_>::
GridUniformQuadCache(shared_ptr<const GridType> grid,
                     const ValueFlags flag,
                     const Quadrature<dim> &quad)
:
                     grid_(grid),
                     flags_(flag),
                     lengths_(grid->get_element_lengths()),
                     quad_(quad)
                     {}



template <int dim_>
void
GridUniformQuadCache<dim_>::
init_element_cache(ElementIterator &elem)
{
    auto acc = elem.get_accessor().elem_values_;
    const auto n_points_direction = quad_.get_num_points_direction();
    const Size n_points = n_points_direction.flat_size();

    if (flags_.fill_points())
    {
        acc.unit_points_ = quad_.get_points();
        flags_.set_points_filled(true);
    }

    if (flags_.fill_w_measures())
    {
        if (acc.w_measure_.size() != n_points)
            acc.w_measure_.resize(n_points);

        acc.unit_weights_ = quad_.get_weights().get_flat_tensor_product();
    }
    else
    {
        Assert(false, ExcMessage("Should not get here"));
        acc.w_measure_.clear() ;
        acc.unit_weights_.clear() ;
    }
}



template <int dim_>
void
GridUniformQuadCache<dim_>::
fill_element_cache(ElementIterator &elem)
{

}



template <int dim_>
void
GridUniformQuadCache<dim_>::
print_info(LogStream &out) const
{
    out.begin_item("Lengths:");
    lengths_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_uniform_quad_cache.inst>
