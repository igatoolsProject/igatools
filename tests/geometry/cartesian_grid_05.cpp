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
 *  Test for the grid_utils grids union
 *
 *  author: martinelli
 *  date: 2013-07-15
 *  QA: (pauletti 2013-10-25)
 *     -use create for grids
 *     -build_cartesian_grid_union should take shared_ptrs not references
 *     -Should this work for <0>, where is this function used?
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/grid_utils.h>

template <int dim>
void run_test()
{
    typedef CartesianGrid<dim> cartesian_grid_t ;

    const int n_knots_grid_0 = 5 ;

    cartesian_grid_t grid_0(n_knots_grid_0);
    grid_0.print_info(out);


    const int n_knots_grid_1 = n_knots_grid_0 - 1 ;
    cartesian_grid_t grid_1(n_knots_grid_1);
    grid_1.print_info(out);


    shared_ptr<cartesian_grid_t> grid_union = build_cartesian_grid_union(grid_0, grid_1) ;
    grid_union->print_info(out);
}


int main()
{
    out.depth_console(10);
    //run_test<0>();
    run_test<1>();
    run_test<2>();
    run_test<3>();

    return  0;
}
