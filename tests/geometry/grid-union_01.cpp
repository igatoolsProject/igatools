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
/*
 *  Test for the grid union: given two different grids grid_1 and grid_2 on the same domain,
 *  it build the grid as the union of grid_1 and grid_2 (i.e. the grid that contains both)
 *
 *  author: martinelli
 *  date: 2014-05-14
 */

// TODO (pauletti, Aug 6, 2014): put some numbear like _01 to filename
#include "../tests.h"

#include <igatools/geometry/grid_tools.h>
#include <igatools/geometry/cartesian_grid_element.h>

template<int dim>
void do_test()
{
    out << "=============== Dimension " << dim << " -- begin ===============" << endl ;

    auto grid_1 = CartesianGrid<dim>::create(3);
    out << "------------------------------------------" << endl;
    out << "Grid 1" << endl;
    grid_1->print_info(out);
    out << "------------------------------------------" << endl;

    out << endl;

    auto grid_2 = CartesianGrid<dim>::create(4);
    out << "------------------------------------------" << endl;
    out << "Grid 2" << endl;
    grid_2->print_info(out);
    out << "------------------------------------------" << endl;

    out << endl;

    grid_tools::InterGridMap map_elem_grid_union_to_elem_grid_1;
    grid_tools::InterGridMap map_elem_grid_union_to_elem_grid_2;
    auto grid_union = grid_tools::build_cartesian_grid_union
                      (*grid_1, *grid_2,
                       map_elem_grid_union_to_elem_grid_1,map_elem_grid_union_to_elem_grid_2);
    out << "------------------------------------------" << endl;
    out << "Grid union" << endl;
    grid_union->print_info(out);

    out << "map_elem_grid_union_to_elem_grid_1 = ";
    for (auto map : map_elem_grid_union_to_elem_grid_1)
        out << map.second << " ";
    out << endl;
    out << "map_elem_grid_union_to_elem_grid_2 = ";
    for (auto map : map_elem_grid_union_to_elem_grid_2)
        out << map.second << " ";
    out << endl;

    out << "------------------------------------------" << endl;
    out << "=============== Dimension " << dim << " -- end =================" << endl ;

    out << endl << endl;
}

int main()
{

    do_test<0>();
    do_test<1>();
    do_test<2>();
    do_test<3>();

    return 0;
}
