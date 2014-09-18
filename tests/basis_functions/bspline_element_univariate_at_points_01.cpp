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
 *  Test for the Bspline space element iterator
 *  Computes univariate values and derivatives of the basis functions
 *
 *  author: martinelli
 *  date: Sep 17, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>


template<int dim, int range >
void do_test()
{
    OUTSTART

    out << "do_test<" << dim << "," << range << ">" << endl;

    auto knots = CartesianGrid<dim>::create();

    const int degree = 2;
    const int rank =  1;
    typedef BSplineSpace< dim, range, rank > Space_t;
    auto space = Space_t::create(degree, knots);

    const int n_points = 3;
    QGauss< dim > quad_scheme(n_points);

    auto element = space->begin();

    vector<string> names= {"Values","Gradients","Hessians"};
    const int max_deriv_order = names.size()-1;

    for (int deriv_order = 0 ; deriv_order <= max_deriv_order ; ++deriv_order)
    {
        auto univariate_values =
            element->evaluate_univariate_derivatives_at_points(deriv_order,quad_scheme);

        out << names[deriv_order] << " 1D:" << endl;

        int comp = 0;
        for (const auto &values_comp : univariate_values)
        {
            out << "Component[" << comp++ << "] : " << endl;
            out.push("  ");
            int dir = 0;
            for (const auto &values_comp_dir : values_comp)
            {
                out << "Direction[" << dir++ << "]" << endl;
                out.push("  ");
                values_comp_dir.print_info(out);
                out.pop();
            }
            out.pop();
        }
        out << endl;
    }

    OUTEND
}


int main()
{
    out.depth_console(10);

    do_test<1,1>();
    do_test<1,2>();
    do_test<1,3>();

    do_test<2,1>();
    do_test<2,2>();
    do_test<2,3>();

    do_test<3,1>();
    do_test<3,3>();

    return 0;
}
