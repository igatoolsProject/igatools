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
 *  Test for the CartesianGrid ElementIterator using UniformQuadCache
 *
 *  author: pauletti
 *  date: 2014-08-15
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid_uniform_quad_cache.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>


template <int dim>
void run_test()
{

    const int n_knots = 5;
    auto grid = CartesianGrid<dim>::create(n_knots);

    QGauss<dim> q1(2);
    GridElementHandler<dim> cache1(grid, ValueFlags::w_measure, q1);

    QGauss<dim> q2(1);
    GridElementHandler<dim> cache2(grid, ValueFlags::w_measure, q2);



    auto el1 = grid->begin();
    cache1.init_element_cache(el1);
    cache2.init_element_cache(el1);

    auto el2 = grid->begin();
    cache2.init_element_cache(el2);
}


int main()
{

    run_test<1>();
    run_test<2>();
    run_test<3>();

    return  0;
}
