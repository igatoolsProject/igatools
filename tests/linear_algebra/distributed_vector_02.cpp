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
 *  Test for Vector class
 *
 *  author: pauletti
 *  date: Apr 5, 2013
 *
 */

#include "../tests.h"

#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/base/ig_function.h>
#include <igatools/base/function_element.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>

template <LAPack la_pack>
void non_contig_indices()
{
	OUTSTART

    using VectorType = Vector<la_pack>;

    vector<Index> dofs = {1, 3, 5};
    VectorType vec(dofs);
    vec.print_info(out);
    
    const int dim = 1;
    using Space = BSplineSpace<dim>;
    using Function = IgFunction<ReferenceSpace<dim>>;
    auto grid = CartesianGrid<dim>::create(3);
    const int deg = 1;
    auto space = Space::create(deg, grid);
    auto F = Function::create(space, vec);

    auto vec1 = F->get_coefficients();
    vec1.print_info(out);

    auto func = F->clone();
    auto vec2 = std::dynamic_pointer_cast<Function>(func)->get_coefficients();
    vec2.print_info(out);

    OUTEND
}



int main()
{
	non_contig_indices<LAPack::trilinos_epetra>();
    return  0;
}
