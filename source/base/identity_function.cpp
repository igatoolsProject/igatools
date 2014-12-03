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

#include <igatools/base/identity_function.h>
#include <igatools/base/function_element.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim,int space_dim>
IdentityFunction<dim,space_dim>::
IdentityFunction(std::shared_ptr<GridType> grid)
    :
    parent_t::NewFunction(grid)
{}



template<int dim,int space_dim>
auto
IdentityFunction<dim,space_dim>::
fill_cache(ElementAccessor &elem, const int j, const variant_2 &k) -> void
{
    parent_t::fill_cache(elem, j, k);
    fill_cache_impl.j = j;
    fill_cache_impl.function = this;
    fill_cache_impl.elem = &elem;
    fill_cache_impl.flags_ = &(this->flags_);
    boost::apply_visitor(fill_cache_impl, k);
}

IGA_NAMESPACE_CLOSE

#include <igatools/base/identity_function.inst>
