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

#include <igatools/basis_functions/bspline_uniform_quad_cache.h>


using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim_, int range_ , int rank_>
BsplineUniformQuadCache<dim_, range_, rank_>::
BsplineUniformQuadCache(shared_ptr<const GridType> grid,
                     const ValueFlags flag,
                     const Quadrature<dim> &quad)
    :
    grid_(grid),
    flags_(flag),
    lengths_(grid->get_element_lengths()),
    quad_(quad)
{}



template<int dim_, int range_ , int rank_>
void
BsplineUniformQuadCache<dim_, range_, rank_>::
init_element_cache(ElementIterator &elem)
{
    // TODO (pauletti, Aug 14, 2014): create get_cache in accessor
    auto &cache = elem.get_accessor().elem_values_;
    cache.resize(flags_, quad_);
}



template<int dim_, int range_ , int rank_>
void
BsplineUniformQuadCache<dim_, range_, rank_>::
fill_element_cache(ElementIterator &elem)
{
    auto &cache = elem.get_accessor().elem_values_;
    auto meas = lengths_.tensor_product(elem->get_tensor_index());
    cache.fill(meas);
    cache.set_filled(true);
}



template<int dim_, int range_ , int rank_>
void
BsplineUniformQuadCache<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out.begin_item("Lengths:");
    lengths_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_uniform_quad_cache.inst>
