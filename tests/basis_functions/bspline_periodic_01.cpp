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
 *  Test for periodic BSplineSpace constructors
 *
 *  author: pauletti
 *  date: 2014-10-23
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/ig_function.h>
#include <igatools/base/identity_function.h>
#include <igatools/io/writer.h>



template <int dim>
void plot_basis(const int n_knots, const int deg)
{
    using Space  = BSplineSpace<dim>;
    using RefSpace  = ReferenceSpace<dim>;
    using Coeffs = typename IgFunction<RefSpace>::CoeffType;

    TensorIndex<dim> deg1(deg);
    typename Space::DegreeTable degt(deg1);

    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto space = Space::create(deg, grid, InteriorReg::maximum, true,
                               BasisEndBehaviour::periodic);


    const int n_basis = space->get_num_basis();
    out << "Num basis: " << n_basis << endl;
    space->print_info(out);

#if 0
    Coeffs coeffs(n_basis);

    for (int basis_index = 0; basis_index < space->get_num_basis(); ++basis_index)
    {
        coeffs[basis_index] = 1.;

        const int n_plot_points = 10;
        Writer<dim> output(IdentityFunction<dim>::create(grid), n_plot_points);

        string field_name = "basis " + to_string(basis_index);

        auto basis = IgFunction<Space>::create(space, coeffs);
        output.template add_field<1,1>(basis, field_name);

        string file_name = "bspline_basis-" + to_string(basis_index);
        output.save(file_name);

        coeffs[basis_index] = 0.;
    }
#endif

}


int main()
{
    const int deg = 1;
    //const int deg = 2;
    const int n_knots = 5 + deg;
    plot_basis<1>(n_knots, deg);
    plot_basis<2>(n_knots, deg);
    plot_basis<3>(n_knots, deg);
    return 0;
}
