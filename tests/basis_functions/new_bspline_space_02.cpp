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
 *  Test for developin new BsplineSpace
 *  author: pauletti
 *  date:
 *
 */
#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>



int main()
{
    out.depth_console(10);

    {
        const int dim=1;
        using BSplineSpace = BSplineSpace<dim>;
        using MultiplicityTable = typename BSplineSpace::MultiplicityTable;
        using DegreeTable = typename BSplineSpace::DegreeTable;

        DegreeTable deg{{2}};
        auto grid = CartesianGrid<dim>::create(4);
        auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable ({ {{1,3}} }));
        auto space = BSplineSpace::create(deg, grid, int_mult);

        space->print_info(out);
        out << endl;
    }

    {
        const int dim=2;
        using BSplineSpace = BSplineSpace<dim>;
        using MultiplicityTable = typename BSplineSpace::MultiplicityTable;
        using DegreeTable = typename BSplineSpace::DegreeTable;

        DegreeTable deg{{2,1}};
        auto grid = CartesianGrid<dim>::create({4,3});
        auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable ({ {{1,3},{1}} }));
        auto space = BSplineSpace::create(deg, grid, int_mult);

        space->print_info(out);
        out << endl;
    }

    {
        const int dim=2;
        const int range=2;
        using BSplineSpace = BSplineSpace<dim,range>;
        using MultiplicityTable = typename BSplineSpace::MultiplicityTable;
        using DegreeTable = typename BSplineSpace::DegreeTable;

        DegreeTable deg{{2,1},{1,3}};
        auto grid = CartesianGrid<dim>::create({4,3});
        auto int_mult = shared_ptr<MultiplicityTable>(
                new MultiplicityTable ({ {{1,3},{1}}, {{1,1},{2}}}));
        auto space = BSplineSpace::create(deg, grid, int_mult);

        space->print_info(out);
        out << endl;
    }


    // Maximum regularity constructor
    {
        const int dim=2;
        using BSplineSpace = BSplineSpace<dim>;
        using DegreeTable = typename BSplineSpace::DegreeTable;

        DegreeTable deg{{2,1}};
        auto grid = CartesianGrid<dim>::create({4,3});
        auto space = BSplineSpace::create(deg, grid);

        space->print_info(out);
        out << endl;
    }

    // Maximum regularity constructor
    {
        const int dim=2;
        using BSplineSpace = BSplineSpace<dim>;

        int deg=3;
        auto grid = CartesianGrid<dim>::create({4,3});
        auto space = BSplineSpace::create(deg, grid);

        space->print_info(out);
        out << endl;
    }

    return 0;
}
