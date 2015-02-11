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
 *  Test for BSplineSpace constructors
 *
 *  author: antolin
 *  date: 2013-12-16
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>

template<int dim>
void def_const(const int &degree_0 = 0)
{

    const int n_knots = 2 ;
    auto grid = CartesianGrid<dim>::create(n_knots) ;
    TensorIndex<dim> degree;
    for (int i = 0; i < dim; ++i)
    {
        degree[i] = degree_0 + i ;
    }
    auto space = BSplineSpace<dim>::create(degree, grid) ;

    out << "Initial degree = " << degree_0 << std::endl;
    space->print_info(out) ;
    out << endl;

}

int main()
{
    def_const<0>(0);
    def_const<1>(1);
    def_const<2>(2);
    def_const<3>(3);

    return 0;
}
