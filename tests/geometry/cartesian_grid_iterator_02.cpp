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
 *  Test for the CartesianGrid element iterator init_value
 *
 *  author: pauletti
 *  date: apr 4, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>


template <int dim>
void run_test()
{

    const int n_knots = 5;
    auto grid = CartesianGrid<dim>::create(n_knots);

    auto el1 = grid->begin();
    el1->init_cache(ValueFlags::w_measure, QGauss<dim>(2));

    el1->init_cache(ValueFlags::w_measure, QGauss<dim>(1));

    auto el2 = grid->begin();
    el2->init_cache(ValueFlags::w_measure, QGauss<dim>(1));
}

int main()
{
    out.depth_console(10);
    run_test<1>();
    run_test<2>();
    run_test<3>();

    return  0;
}
