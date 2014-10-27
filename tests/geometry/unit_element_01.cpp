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
#include <igatools/base/quadrature.h>

template<std::size_t... I>
auto tuple_of_quads(std::index_sequence<I...>)
-> decltype(std::make_tuple(Quadrature<I>() ...))
{
    return std::make_tuple(Quadrature<I>() ...);
}

template<int dim>
using QuadList = decltype(tuple_of_quads(std::make_index_sequence<dim+1>()));



template<class ValuesCache, int dim, std::size_t... I>
auto tuple_of_caches(std::index_sequence<I...>, const Quadrature<dim> &q, const ValuesCache &)
-> decltype(std::make_tuple(std::array<ValuesCache,
                            UnitElement<dim>::template num_elem<dim-I>()>() ...))
{
    return std::make_tuple(std::array<ValuesCache,
                           UnitElement<dim>::template num_elem<dim-I>()>() ...);
}


template<int dim, int n_sub_elem>
using CacheList = decltype(tuple_of_caches(std::make_index_sequence<n_sub_elem+1>(),
                                           Quadrature<dim>(),
                                           3));

int main()
{


    CacheList<3,1> list_values;
    auto &cache0 = std::get<0>(list_values);
    auto &cache1 = std::get<1>(list_values);

    for (auto &c : cache0)
        out << c << endl;

    for (auto &c : cache1)
        out << c << endl;

    UnitElement<2>::num_elem<1>();
    UnitElement<2>::get_elem<1>(0);

    const int dim=3;
    out.depth_console(20);
    QuadList<dim>  list_of_quad;

    auto &quad0 = std::get<0>(list_of_quad);
    quad0.print_info(out);
    out << endl;

    auto &quad1 = std::get<1>(list_of_quad);
    quad1.print_info(out);
    out << endl;


    auto &quad2 = std::get<2>(list_of_quad);
    quad2.print_info(out);
    out << endl;
    return 0;
}
