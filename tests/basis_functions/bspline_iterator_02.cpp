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
 *  Test for the BSplineSpace element iterator local to global
 *  author: pauletti
 *  date: Aug 21, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>



template <int dim, int range=1, int rank=1>
void elem_dofs(const int n_knots = 4, const int deg=1)
{
    OUTSTART

    using Space = NewBSplineSpace<dim, range, rank>;
    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto space = Space::create(deg, grid);

    auto elem = space->begin();
    auto end = space->end();

    for (; elem != end; ++elem)
    {
        out << "Element index: " << elem->get_flat_index() << endl;
        out << "Global dofs: ";
        elem->get_local_to_global().print_info(out);
        out << endl;
    }

    OUTEND
}


int main()
{
    out.depth_console(10);

    elem_dofs< 1, 1, 1 >();
    elem_dofs< 1, 2, 1 >();
    elem_dofs< 1, 3, 1 >();
    elem_dofs< 2, 1, 1 >();
    elem_dofs< 2, 2, 1 >();
    elem_dofs< 2, 3, 1 >();
    elem_dofs< 3, 1, 1 >();
    elem_dofs< 3, 3, 1 >();

    return  0;
}
