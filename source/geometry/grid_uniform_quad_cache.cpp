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
using std::array;

IGA_NAMESPACE_OPEN

template <int dim_>
const array<Size, UnitElement<dim_>::faces_per_element>
GridUniformQuadCache<dim_>::faces  = UnitElement<dim_>::faces;

template <int dim_>
GridUniformQuadCache<dim_>::
GridUniformQuadCache(shared_ptr<const GridType> grid,
                     const ValueFlags flag,
                     const Quadrature<dim> &quad)
    :
    grid_(grid),
    flags_(flag),
    face_flags_(flag),
    lengths_(grid->get_element_lengths()),
    quad_(quad)
{}



template <int dim_>
void
GridUniformQuadCache<dim_>::
init_element_cache(ElementAccessor &elem)
{
    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    auto &elem_cache = cache->elem_values_;
    elem_cache.resize(flags_, quad_);

    auto &face_cache = cache->face_values_;
    for (auto &f: faces)
        face_cache[f].resize(face_flags_, quad_, f);
}



template <int dim_>
void
GridUniformQuadCache<dim_>::
init_element_cache(ElementIterator &elem)
{
    init_element_cache(elem.get_accessor());
}


template <int dim_>
void
GridUniformQuadCache<dim_>::
fill_element_cache(ElementAccessor &elem)
{
    const auto &index = elem.get_tensor_index();
    auto meas = lengths_.tensor_product(index);

    Assert(elem.local_cache_ != nullptr,ExcNullPtr());
    auto &cache = elem.local_cache_->elem_values_;
    cache.fill(meas);
    cache.set_filled(true);
}



template <int dim_>
void
GridUniformQuadCache<dim_>::
fill_element_cache(ElementIterator &elem)
{
    fill_element_cache(elem.get_accessor());
}



template <int dim_>
void
GridUniformQuadCache<dim_>::
fill_face_cache(ElementIterator &elem, const int face)
{
    const auto &index = elem->get_tensor_index();

    auto &elem_accessor = elem.get_accessor();
    Assert(elem_accessor.local_cache_ != nullptr,ExcNullPtr());
    auto &f_cache = elem_accessor.local_cache_->face_values_[face];
    auto meas = lengths_.sub_tensor_product(index, UnitElement<dim_>::face_active_directions[face]);
    f_cache.fill(meas);
    f_cache.set_filled(true);
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
