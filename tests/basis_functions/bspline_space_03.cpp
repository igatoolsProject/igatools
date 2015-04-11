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
 *  Test for building a matrix on a space of an igfunction
 *
 *  author: pauletti
 *  date: 2015-03-17
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/ig_function.h>
#include <igatools/linear_algebra/epetra_matrix.h>


template<int dim, int codim=0>
void using_const_space(shared_ptr<IgFunction<ReferenceSpace<dim>>> fun)
{
    OUTSTART

    auto space = fun->get_ig_space();
    Epetra_SerialComm comm;
    auto map = EpetraTools::create_map(space, "active", comm);
    auto graph = EpetraTools::create_graph(space, "active", space, "active", map, map);
    auto matrix = EpetraTools::create_matrix(graph);

    OUTEND
}

template<int dim, int codim=0>
void using_const_space(shared_ptr<Function<dim> fun)
{
    OUTSTART
    auto grid = fun->get_hacked_grid();
    auto zero = Function::create(grid,
                                IdentityFunction<dim>::create(grid), val);
    OUTEND
}

int main()
{
    const int dim = 2;
    using Space = BSplineSpace<dim>;

    auto grid = CartesianGrid<dim>::create(5);
    auto space = Space::create(1, grid);

    auto coeff = EpetraTools::create_vector(space, "active");

    auto fun = IgFunction<ReferenceSpace<dim>>::create(space, *coeff);

    using_const_space<2>(fun);

    return 0;
}
