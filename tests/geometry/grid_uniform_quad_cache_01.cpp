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
 *  Test for the CartesianGrid UniformQuadCache
 *
 *  author: pauletti
 *  date: Aug 21, 2014
 *
 */

#include "../tests.h"

#include <igatools/geometry/grid_uniform_quad_cache.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>


template <int dim>
void run_test()
{
    out << "dim = " << dim << endl;
    const int n_knots = 5;
    auto grid = CartesianGrid<dim>::create(n_knots);
    auto flag = ValueFlags::measure|ValueFlags::w_measure;
    auto quad = QGauss<dim>(2);
    GridElementHandler<dim> cache(grid, flag, quad);
    cache.print_info(out);
    out << endl;

    auto elem = grid->begin();
    cache.init_element_cache(elem);

    auto end = grid->end();
    for (; elem != end; ++elem)
    {
        cache.fill_element_cache(elem);
        out << "Measure: " << elem->get_measure() << endl;
    }
}


int main()
{
    out.depth_console(10);
    run_test<1>();
    run_test<2>();
    run_test<3>();

    return  0;
}
