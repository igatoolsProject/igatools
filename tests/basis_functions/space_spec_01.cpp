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
 *  Test for the SplineSpace class
 *  author: pauletti
 *  date:
 *
 */

#include "../tests.h"
#include <igatools/basis_functions/spline_space.h>

using std::shared_ptr;

void test_1d()
{
    const int dim=1;
    using SplineSpace = SplineSpace<dim>;
    using MultiplicityTable = typename SplineSpace::MultiplicityTable;

    auto grid = CartesianGrid<dim>::create(4);
    typename SplineSpace::DegreeTable deg {{2}};
    auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable({ {{1,3}} }));
    SplineSpace sp_spec(deg, grid, int_mult);
    sp_spec.print_info(out);

    CartesianProductArray<Real,2> bn_x {{-0.5, 0, 0}, {1.1, 1.2, 1.3}};
    typename SplineSpace::BoundaryKnotsTable bdry_knots { {bn_x} };
    auto rep_knots = sp_spec.compute_knots_with_repetition(bdry_knots);
    out << "Boundary knots:\n";
    for (const auto &v : bdry_knots)
        for (const auto &w : v)
            w.print_info(out);
    out << "Repeated knots:\n";
    for (const auto &v : rep_knots)
        v.print_info(out);
}



void test_2d()
{
    const int dim=2;
    using SplineSpace = SplineSpace<dim>;
    using MultiplicityTable = typename SplineSpace::MultiplicityTable;
    auto grid = CartesianGrid<dim>::create({3,5});
    typename SplineSpace::DegreeTable deg {{1,3}};

    auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable({ {{1}, {1,3,1}} }));

    SplineSpace sp_spec(deg, grid, int_mult);
    sp_spec.print_info(out);

    iga::CartesianProductArray<double, 2> bk_x {{-0.5, 0}, {1.2, 1.3}};
    iga::CartesianProductArray<double, 2> bk_y {{-0.6,0,0,0}, {1,1.1,1.6, 1.6}};
    typename SplineSpace::BoundaryKnotsTable bdry_knots { {bk_x, bk_y} };
    auto rep_knots = sp_spec.compute_knots_with_repetition(bdry_knots);
    out << "Boundary knots:\n";
    for (const auto &v : bdry_knots)
        for (const auto &w : v)
            w.print_info(out);
    out << "Repeated knots:\n";
    for (const auto &v : rep_knots)
        v.print_info(out);
}


void test_3d()
{
    const int dim=3;
    using SplineSpace = SplineSpace<dim>;
    using MultiplicityTable = typename SplineSpace::MultiplicityTable;
    auto grid = CartesianGrid<dim>::create({3,4,5});
    typename SplineSpace::DegreeTable deg {{1,3,0}};
    auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable({ {{1}, {1,3}, {1,1,1}} }));




    SplineSpace sp_spec(deg, grid, int_mult);
    sp_spec.print_info(out);

    iga::CartesianProductArray<double, 2> bk_x {{-0.5, 0}, {1.2, 1.3}};
    iga::CartesianProductArray<double, 2> bk_y {{-0.6,0,0,0}, {1,1,1.6, 1.6}};
    iga::CartesianProductArray<double, 2> bk_z {{-0.6}, {1.6}};
    typename SplineSpace::BoundaryKnotsTable bdry_knots { {bk_x, bk_y, bk_z} };

    auto rep_knots = sp_spec.compute_knots_with_repetition(bdry_knots);
    out << "Boundary knots:\n";
    for (const auto &v : bdry_knots)
        for (const auto &w : v)
            w.print_info(out);
    out << "Repeated knots:\n";
    for (const auto &v : rep_knots)
        v.print_info(out);
}


void test_2d_2()
{
    const int dim=2;
    const int range=2;
    using SplineSpace = SplineSpace<dim, range, 1>;
    using MultiplicityTable = typename SplineSpace::MultiplicityTable;
    auto grid = CartesianGrid<dim>::create({3,4});
    typename SplineSpace::DegreeTable deg {{1,3},{3,1}};

    auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable({ {{1}, {1,3}},{{1}, {1,1}}}));

    SplineSpace sp_spec(deg, grid, int_mult);
    sp_spec.print_info(out);

    iga::CartesianProductArray<double, 2> bk_x {{-0.5, 0}, {1.2, 1.3}};
    iga::CartesianProductArray<double, 2> bk_y {{-0.6,0,0,0}, {1,1,1.6, 1.6}};

    typename SplineSpace::BoundaryKnotsTable bdry_knots { {bk_x, bk_y}, {bk_y, bk_x} };
    auto rep_knots = sp_spec.compute_knots_with_repetition(bdry_knots);
    out << "Boundary knots:\n";
    for (const auto &v : bdry_knots)
        for (const auto &w : v)
            w.print_info(out);
    out << "Repeated knots:\n";
    for (const auto &v : rep_knots)
        v.print_info(out);
}


int main()
{
    out.depth_console(10);

    test_1d();
    test_2d();
    test_3d();

    test_2d_2();

    return 0;
}
