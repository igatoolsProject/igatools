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

#include <igatools/basis_functions/new_bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
// [new include]
#include <igatools/base/ig_function.h>

// TODO (pauletti, Nov 26, 2014): this is not correct, fix the writer
#include <igatools/base/identity_function.h>
// [new include]
#include <igatools/io/writer.h>
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

// [plot function]
template <int dim>
void plot_basis(const int deg)
{
    using Space  = BSplineSpace<dim>;
    using Coeffs = typename IgFunction<Space>::CoeffType;

    const int n_knots = deg + 2;
    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto space = Space::create(deg, grid);
    // [plot function]

    // [init vec]
    const int n_basis = space->get_num_basis();
    Coeffs coeffs(n_basis);
    // [init vec]

    // [tensor to flat]
    TensorIndex<dim> basis_t_index(deg);
    auto basis_index = space->get_global_dof_id(basis_t_index, 0);
    coeffs[basis_index] = 1.;
    // [tensor to flat]

    // [print vector]
    out << "Coefficient vector of: " << basis_index << "-th basis" << endl;
    coeffs.print_info(out);
    out << endl;
    // [print vector]

    // [plot basis]
    out << "Saving basis plot" << endl;
    const int n_plot_points = 5;
    Writer<dim> output(IdentityFunction<dim>::create(grid), n_plot_points);

    string field_name = "basis " + to_string(basis_index);

    auto basis = IgFunction<Space>::create(space, coeffs);
    output.template add_field<1,1>(basis, field_name);

    string file_name = "bspline_basis-" + to_string(dim) + "d";
    output.save(file_name);
    // [plot basis]

}


int main()
{
    const int deg = 2;

    plot_basis<1>(deg);
    plot_basis<2>(deg);
    plot_basis<3>(deg);

    return 0;
}



