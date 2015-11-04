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

#include <igatools/basis_functions/bspline.h>
#include <igatools/io/writer.h>
// [new include]
#include <igatools/functions/ig_grid_function.h>
// [new include]

using namespace iga;
using std::endl;
using std::to_string;

LogStream out;

// [plot_function]
template <int dim>
void plot_basis(const int deg)
{
  const int n_knots = deg + 2;
  const auto grid = Grid<dim>::const_create(n_knots);
// [plot_function]

// [create_basis]
  const auto space = SplineSpace<dim>::const_create(deg, grid);
  const auto basis = BSpline<dim>::const_create(space);
// [create_basis]

  // [init_vec]
  IgCoefficients coeffs(basis->get_global_dofs());
  // [init_vec]

  //[tensor_to_flat]
  TensorIndex<dim> basis_t_index(deg);
  const auto j = basis->get_global_dof_id(basis_t_index, 0);
  coeffs[j] = 1.0;
  // [tensor_to_flat]

  // [print_vector]
  out << "IgCoefficient for: " << j << "-th basis" << endl;
  coeffs.print_info(out);
  out << endl;
  // [print_vector]

  // [basis_to_plot]
  auto central_basis = IgGridFunction<dim,1>::const_create(basis,coeffs);
  // [basis_to_plot]

  // [plot_basis]
  out << "Saving basis plot" << endl;
  const int n_plot_points = 10;
  Writer<dim> output(grid,n_plot_points);

  string func_name = "basis " + to_string(j);
  output.template add_field(*central_basis, func_name);

  string file_name = "bspline_basis-" + to_string(j) + "_" + to_string(dim) + "d";
  output.save(file_name);
  // [plot_basis]
}


int main()
{
  const int deg = 2;

  plot_basis<1>(deg);
  plot_basis<2>(deg);
  plot_basis<3>(deg);

  return 0;
}



