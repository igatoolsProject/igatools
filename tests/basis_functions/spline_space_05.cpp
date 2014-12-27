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
 *  Test for the SplineSpace class homogeneous range
 *  author: pauletti
 *  date:
 *
 */

#include "../tests.h"
#include <igatools/basis_functions/spline_space.h>

template<int dim, int range, int rank>
void test(const int deg1)
{
    using SplineSpace = SplineSpace<dim, range, rank>;

    auto grid = CartesianGrid<dim>::create(4);
    typename SplineSpace::Degrees deg2(deg1);
    typename SplineSpace::DegreeTable deg(deg2);

    deg.get_active_components_id().print_info(out);
    out << endl;
    deg.get_inactive_components_id().print_info(out);
    out << endl;

    auto int_mult = SplineSpace::multiplicity_regularity(InteriorReg::maximum,
                                                         deg, grid->get_num_intervals());
    SplineSpace sp_spec(deg, grid, int_mult);
    sp_spec.print_info(out);
}



int main()
{
    out.depth_console(10);

    test<1, 1, 1>(1);
    test<1, 2, 1>(1);
    test<2, 2, 1>(3);
    test<2, 2, 1>(2);

    return 0;
}
