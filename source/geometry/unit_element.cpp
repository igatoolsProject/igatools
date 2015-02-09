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

#include <igatools/geometry/unit_element.h>
#include <igatools/base/array_utils.h>

IGA_NAMESPACE_OPEN

template <int dim>
std::array<Size, dim+1>
fill_skeleton_size()
{
    std::array<Size,dim+1> res;
    for (int k=0; k<dim+1; ++k)
        res[k] = skel_size(dim, k);
    return res;
}


template <int dim>
const std::array<Size, dim>
UnitElement<dim>::active_directions = sequence<dim>();


template <int dim>
const std::array<Size, dim + 1>
UnitElement<dim>::sub_elements_size = fill_skeleton_size<dim>();



template <int dim>
const decltype(tuple_of_elements<dim>(std::make_index_sequence<dim+1>()))
UnitElement<dim>::all_elems = construct_cube_elements<dim>();

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/unit_element.inst>
