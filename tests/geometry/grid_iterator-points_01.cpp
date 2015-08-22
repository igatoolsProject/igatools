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
 *  Test for the CartesianGrid ElementIterator get_points()
 *
 *  author: pauletti
 *  date: 2014-08-15
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/grid_cache_handler.h>
#include <igatools/geometry/grid_element.h>


template <int dim>
void elem_points(const int n_knots = 5)
{
    OUTSTART

    using Grid = CartesianGrid<dim>;
    using ElementHandler = typename Grid::ElementHandler;

    auto grid = Grid::create(n_knots);

    auto flag = ValueFlags::point;
    QGauss<dim> quad(2);
    ElementHandler cache(grid);
    cache.template reset<dim>(flag, quad);
    auto elem = grid->begin();
    cache.init_element_cache(elem);

    for (; elem != grid->end(); ++elem)
    {
        elem->print_info(out);

        cache.fill_element_cache(elem);
        out.begin_item("Points:");
        elem->get_points().print_info(out);
        out.end_item();
    }

    OUTEND
}


int main()
{
    elem_points<1>();
    elem_points<2>();
    elem_points<3>();

    return  0;
}
