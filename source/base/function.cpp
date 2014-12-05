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

#include <igatools/base/function.h>
#include <igatools/base/function_element.h>

IGA_NAMESPACE_OPEN

template<int dim, int codim, int range, int rank >
Function<dim, codim, range, rank >::
Function(std::shared_ptr<GridType> grid)
    :
    GridElementHandler<dim>(grid)
{}



template<int dim, int codim, int range, int rank>
auto
Function<dim, codim, range, rank >::
get_cache(ElementAccessor &elem)
-> std::shared_ptr<typename ElementAccessor::CacheType> &
{
    return elem.local_cache_;
}



IGA_NAMESPACE_CLOSE

#include <igatools/base/function.inst>

