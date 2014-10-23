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
 *  Test for BSplineSpace constructors
 *
 *  author: pauletti
 *  date: 2014-10-23
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/new_bspline_space.h>

namespace grid
{
template<int dim>
shared_ptr<CartesianGrid<dim>>
uniform(const int n_knots)
{
    return CartesianGrid<dim>::create(n_knots);
}


};


template<int dim>
void uniform_degree(const int deg, shared_ptr<CartesianGrid<dim>> grid)
{
    OUTSTART
    auto space = NewBSplineSpace<dim>::create(deg, grid);
    space->print_info(out);
    out << endl;
    OUTEND
}


template<int dim>
void direction_degree(const TensorIndex<dim> &deg,
                      shared_ptr<CartesianGrid<dim>> grid)
{
    OUTSTART
    auto space = NewBSplineSpace<dim>::create(deg, grid);
    space->print_info(out);
    out << endl;
    OUTEND
}


int main()
{
    const int deg = 1;
    const int n_knots = 2;
    uniform_degree<0>(deg, grid::uniform<0>(n_knots));
    uniform_degree<1>(deg, grid::uniform<1>(n_knots));
    uniform_degree<2>(deg, grid::uniform<2>(n_knots));
    uniform_degree<3>(deg, grid::uniform<3>(n_knots));

    TensorIndex<0> deg0;
    direction_degree<0>(deg0, grid::uniform<0>(n_knots));

    TensorIndex<1> deg1 = {1};
    direction_degree<1>(deg1, grid::uniform<1>(n_knots));

    TensorIndex<2> deg2 = {2,3};
    direction_degree<2>(deg2, grid::uniform<2>(n_knots));

    TensorIndex<3> deg3 = {3,4,5};
    direction_degree<3>(deg3, grid::uniform<3>(n_knots));

    return 0;
}
