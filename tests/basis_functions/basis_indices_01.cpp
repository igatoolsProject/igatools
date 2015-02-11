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
 *  Test for basis index
 *  author: pauletti
 *  date:
 *
 */

#include "../tests.h"
#include <igatools/basis_functions/dof_distribution.h>

int main()
{
    out.depth_console(10);

    {
        const int dim=1;
        using SplineSpace = SplineSpace<dim>;
        using MultiplicityTable = typename SplineSpace::MultiplicityTable;

        typename SplineSpace::DegreeTable deg {{2}};

        auto grid = CartesianGrid<dim>::create(4);

        auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable({ {{1,3}} }));
        SplineSpace sp_spec(deg, grid, int_mult);

        CartesianProductArray<Real,2> bn_x {{-0.5, 0, 0}, {1.1, 1.2, 1.3}};
        typename SplineSpace::BoundaryKnotsTable bdry_knots { {bn_x} };
        auto rep_knots = sp_spec.compute_knots_with_repetition(bdry_knots);
        auto acum_mult = sp_spec.accumulated_interior_multiplicities();

        auto n_basis = sp_spec.get_num_basis_table();
        auto degree = sp_spec.get_degree();

        DofDistribution<dim> basis_index(grid, acum_mult, n_basis, degree);
        basis_index.print_info(out);
    }

    {
        const int dim=2;
        using SplineSpace = SplineSpace<dim>;

        typename SplineSpace::DegreeTable deg {{1,2}};

        auto grid = CartesianGrid<dim>::create({4,3});

        SplineSpace sp_spec(deg,grid, SplineSpace::InteriorReg::maximum);


        auto rep_knots = sp_spec.compute_knots_with_repetition(sp_spec.get_end_behaviour());
        auto acum_mult = sp_spec.accumulated_interior_multiplicities();

        auto n_basis = sp_spec.get_num_basis_table();
        auto degree = sp_spec.get_degree();

        DofDistribution<dim> basis_index(grid, acum_mult, n_basis, degree);
        basis_index.print_info(out);
    }

    return 0;
}
