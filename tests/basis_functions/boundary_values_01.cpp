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
 *  Test for the boundary l2 projection function.
 *  On a PhysicalSpace
 *
 *  author: pauletti
 *  date: 2014-11-14
 *
 */

#include "../tests.h"
#include "common_functions.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/linear_algebra/dof_tools.h>

#include <igatools/base/identity_function.h>



template<int dim , int codim, int range ,int rank, LAPack la_pack>
void do_test(const int p, const int num_knots = 10)
{
    const int sub_dim = dim - 1;
    using RefSpace = ReferenceSpace<dim,range,rank>;
    using BspSpace = BSplineSpace<dim,range,rank>;
    using Space = PhysicalSpace<RefSpace, codim, Transformation::h_grad>;


    auto grid = CartesianGrid<dim>::create(num_knots);
    auto ref_space = BspSpace::create(p, grid) ;
    auto map = IdentityFunction<dim>::create(grid);
    auto space = Space::create(ref_space, map);

    auto f = BoundaryFunction<dim>::create(grid);


    const int n_qpoints = 4;
    QGauss<sub_dim> quad(n_qpoints);

    const boundary_id dirichlet = 1;
    grid->set_boundary_id(0, dirichlet);
    std::set<boundary_id> bdry_ids;
    bdry_ids.insert(dirichlet);



    std::map<Index,Real> boundary_values;
    space_tools::project_boundary_values<Space,la_pack>(
        f, space, quad, bdry_ids,
        boundary_values);

    out << "basis index \t value" << endl;
    for (auto entry : boundary_values)
        out << entry.first << "\t" << entry.second << endl;

}



int main()
{
#if defined(USE_TRILINOS)
    const auto la_pack = LAPack::trilinos_tpetra;
#elif defined(USE_PETSC)
    const auto la_pack = LAPack::petsc;
#endif
    out.depth_console(20);

    do_test<2, 0, 1, 1, la_pack>(3);

    //do_test<3,3,1,1>(3);

    //do_test<2,3,1,0>(3);

    return 0;
}

