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
#if 0
#include <igatools/geometry/mapping_uniform_quad_cache.h>

using std::shared_ptr;
using std::array;

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
MappingUniformQuadCache<dim_, codim_>::
MappingUniformQuadCache(std::shared_ptr<const Map> map,
                        const ValueFlags flag,
                        const Quadrature<dim> &quad)
    :
    base_t(map->get_grid(), flag, quad),
    flags_(flag),
    quad_(quad)
{}



template<int dim_, int codim_>
void
MappingUniformQuadCache<dim_, codim_>::
init_element_cache(ElementAccessor &elem)
{
    base_t::init_element_cache(elem);
    elem.init_cache(flags_, quad_);
}



template<int dim_, int codim_>
void
MappingUniformQuadCache<dim_, codim_>::
fill_element_cache(ElementAccessor &elem)
{
    base_t::fill_element_cache(elem);
    elem.fill_cache();
}



template<int dim_, int codim_>
void
MappingUniformQuadCache<dim_, codim_>::
init_element_cache(ElementIterator &elem)
{
    init_element_cache(elem.get_accessor());
}



template<int dim_, int codim_>
void
MappingUniformQuadCache<dim_, codim_>::
fill_element_cache(ElementIterator &elem)
{
    fill_element_cache(elem.get_accessor());
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/mapping_uniform_quad_cache.inst>
#endif
