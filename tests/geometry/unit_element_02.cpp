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

/*
 *  Test of unit element
 *  author: pauletti
 *  date: 2014-08-18
 *
 */

#include "../tests.h"

#include <igatools/geometry/unit_element.h>



template <int dim>
void skeleton()
{
    OUTSTART
    for (int k = 0; k<=dim; ++k)
    	out << UnitElement<dim>::skeleton_size[k] << endl;
    OUTEND
}


template <int dim, int k>
EnableIf< (dim==0) || (k<0),
std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)>>
fill_cube_elements()
{
    std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)> res;
    return res;
}

template <int dim, int k>
EnableIf< (dim==k) && (k>0),
std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)>>
fill_cube_elements()
{
    std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)> res;
    res[0].active_directions = sequence<k>();
    return res;
}


template <int dim, int k>
EnableIf< (dim>k) && (k>=0),
std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)>>
fill_cube_elements()
{
    std::array<typename UnitElement<dim>::template Skeleton<k>, skel_size(dim, k)> elements;

    auto sub_elems_1 = fill_cube_elements<dim-1, k>();
    auto sub_elems_0 = fill_cube_elements<dim-1, k-1>();

    auto elem = elements.begin();

    for (auto &sub_elem_0 : sub_elems_0)
    {
        auto &sub_dirs_0 = sub_elem_0.constant_directions;
        auto &dirs       = elem->constant_directions;
        std::copy(sub_dirs_0.begin(), sub_dirs_0.end(), dirs.begin());
        elem->constant_values = sub_elem_0.constant_values;
        ++elem;
    }

    for (int j = 0; j<2; ++j)
    for (auto &sub_elem_1 : sub_elems_1)
    {
        auto &sub_dirs_1 = sub_elem_1.constant_directions;
        auto &sub_values_1 = sub_elem_1.constant_values;
        {
            auto &dirs       = elem->constant_directions;
            auto &values       = elem->constant_values;
            std::copy(sub_dirs_1.begin(), sub_dirs_1.end(), dirs.begin());
            dirs[dim - k -1] = dim-1;
            std::copy(sub_values_1.begin(), sub_values_1.end(), values.begin());
            values[dim - k -1] = j;
            ++elem;
        }
    }

    for (auto &elem : elements)
    {
        auto all = sequence<dim>();
        std::set_difference(all.begin(), all.end(),
                            elem.constant_directions.begin(),
                            elem.constant_directions.end(),
                            elem.active_directions.begin());
    }

    return elements;
}


template<int dim, int sub_dim>
void cube_elements()
{
    OUTSTART
	auto elements = fill_cube_elements<dim, sub_dim>();
    const auto size = elements.size();
    out << "Number of elements: " << size << endl;
	for (auto i=0; i<size; ++i)
	{
	    out.begin_item("Element: " + std::to_string(i));
	    auto &face = elements[i];
	    const auto n_dir = face.constant_directions.size();
        out << "constant directions" << endl;
        for (int j=0; j<n_dir; ++j)
        {
            out << "x["<<face.constant_directions[j]<< "]";
            out << " = " << face.constant_values[j] << endl;
        }
        out << "Active directions" << endl;
        for (auto &dir : face.active_directions)
            out << "x[" << dir << "]" << endl;
        out.end_item();
	}
	OUTEND
}


template<std::size_t... I>
auto tuple_of_elements(std::index_sequence<I...>)
    -> decltype(std::make_tuple(std::array<int, I>() ...))
{
    return std::make_tuple(std::array<int, I>() ...);
}

template<std::size_t N, typename Indices = std::make_index_sequence<N>>
auto cube_elements()
    -> decltype(tuple_of_elements(Indices()))
{
    return tuple_of_elements(Indices());
}

int main()
{
    out.depth_console(20);

    auto all_elems = cube_elements<4>();

    auto faces =  std::get<2>(all_elems);
    for (auto &dir : faces)
    	out << dir << endl;

    return 0;
}
