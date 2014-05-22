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
 *  Test for bernstein extraction class
 *  author: pauletti
 *  date:
 *
 */

#include "../tests.h"
#include <igatools/basis_functions/bernstein_extraction.h>


int main()
{
    out.depth_console(10);

    {
        const int dim=1;
        using SpaceSpec = SpaceSpec<dim>;
        using MultiplicityTable = typename SpaceSpec::MultiplicityTable;

        typename SpaceSpec::DegreeTable deg{{2}};

        auto grid = CartesianGrid<dim>::create(4);

        auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable ({ {{1,3}} }));
        SpaceSpec sp_spec(grid, int_mult, deg);

        CartesianProductArray<Real,2> bn_x{{-0.5, 0, 0}, {1.1, 1.2, 1.3}};
        typename SpaceSpec::BoundaryKnotsTable bdry_knots{ {bn_x} };
        auto rep_knots = sp_spec.compute_knots_with_repetition(bdry_knots);
        auto acum_mult = sp_spec.compute_elements_index_space_mark();


        BernsteinExtraction<dim> operators(grid, rep_knots, acum_mult, deg);
        operators.print_info(out);
    }

    {
        const int dim=1;
        using SpaceSpec = SpaceSpec<dim>;
        using MultiplicityTable = typename SpaceSpec::MultiplicityTable;

        typename SpaceSpec::DegreeTable deg{{3}};

        CartesianProductArray<Real,dim> knots({{0,1,2,3,4}});
        auto grid = CartesianGrid<dim>::create(knots);

        SpaceSpec sp_spec(grid, SpaceSpec::InteriorReg::maximum, deg);


        auto rep_knots = sp_spec.compute_knots_with_repetition(SpaceSpec::EndBehaviour::interpolatory);
        auto acum_mult = sp_spec.compute_elements_index_space_mark();


        BernsteinExtraction<dim> operators(grid, rep_knots, acum_mult, deg);
        operators.print_info(out);
    }

    return 0;
}
