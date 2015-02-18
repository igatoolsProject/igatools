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
 *  Test for active grid elements
 *
 *  author: pauletti
 *  date: Aug 29, 2014
 *
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element.h>

template<int dim>
void
test()
{
    using Grid = CartesianGrid<dim>;
    const int n_knots = 4;
    auto grid = Grid::create(n_knots);
    grid->print_info(out);

    for (auto elem : *grid)
        out << elem.get_flat_index() << endl;

    std::string active_tag = "active";
    grid->add_elements_property(active_tag);

    if (dim > 0)
    {
        std::set<Index> id_all_elems = grid->get_elements_id_same_property(Grid::elems_property_none);


        for (const auto elem_id : id_all_elems)
        {
            if (elem_id % 2 == 0)
                grid->set_element_property_status(active_tag,elem_id,false);
        }
    }

    grid->print_info(out);
    for (auto elem : *grid)
        out << elem.get_flat_index() << endl;
}




int main()
{

    out.depth_console(10);

    test<0>();
    test<1>();
    test<2>();
    test<3>();

    return 0;
}
