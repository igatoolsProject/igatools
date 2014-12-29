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
 *  Test for bernstein extraction class for periodic spaces
 *  author: pauletti
 *  date: 2014-12-22
 *
 */

#include "../tests.h"
#include <igatools/basis_functions/bernstein_extraction.h>

// TODO (pauletti, Dec 26, 2014): this test needs to be update to current standards

template <int dim>
void
test(const int deg = 1)
{
    using SplineSpace = SplineSpace<dim>;
    using MultiplicityTable = typename SplineSpace::MultiplicityTable;

    typename SplineSpace::DegreeTable degt {{deg}};


    CartesianProductArray<Real,dim> knots({{0,1,2,3,4}});
    auto grid = CartesianGrid<dim>::create(knots);

    SplineSpace sp_spec(degt, grid, InteriorReg::maximum, EndBehaviour::periodic);


    auto rep_knots = sp_spec.compute_knots_with_repetition(sp_spec.get_end_behaviour());
    auto acum_mult = sp_spec.accumulated_interior_multiplicities();

    rep_knots.print_info(out);
    out << endl;
    acum_mult.print_info(out);
    out << endl;
    BernsteinExtraction<dim> operators(grid, rep_knots, acum_mult, degt);
    operators.print_info(out);
}


int main()
{
    out.depth_console(10);
    test<1>();


//    {
//        const int dim=1;
//        using SplineSpace = SplineSpace<dim>;
//        using MultiplicityTable = typename SplineSpace::MultiplicityTable;
//
//        typename SplineSpace::DegreeTable deg {{2}};
//
//        auto grid = CartesianGrid<dim>::create(4);
//
//        auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable({ {{1,3}} }));
//        SplineSpace sp_spec(deg, grid, int_mult);
//
//        CartesianProductArray<Real,2> bn_x {{-0.5, 0, 0}, {1.1, 1.2, 1.3}};
//        typename SplineSpace::BoundaryKnotsTable bdry_knots { {bn_x} };
//        auto rep_knots = sp_spec.compute_knots_with_repetition(bdry_knots);
//        auto acum_mult = sp_spec.accumulated_interior_multiplicities();
//
//
//        BernsteinExtraction<dim> operators(grid, rep_knots, acum_mult, deg);
//        operators.print_info(out);
//    }


    return 0;
}
