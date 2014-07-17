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
 * Testing the writer, this the cell data
 * author: pauletti
 * date: Jun 21, 2014
 *
 */

#include "../tests.h"
#include "igatools/io/writer.h"

template<int dim>
void
test()
{
    const int n_knots = 4;
    auto grid = CartesianGrid<dim>::create(n_knots);
    Writer<dim> writer(grid);

    const Size n_iga_elements = writer.get_num_iga_elements();
    const Size n_points_iga_element = writer.get_num_points_per_iga_element();

    vector<vector<vector<Real>>> data_scalar(n_iga_elements,
      vector<vector<Real>>(n_points_iga_element, vector<Real>(1)));

    vector<vector<vector<Real>>> data_vector(n_iga_elements,
      vector<vector<Real>>(n_points_iga_element, vector<Real>(dim)));

    vector<vector<vector<Real>>> data_tensor(n_iga_elements,
      vector<vector<Real>>(n_points_iga_element, vector<Real>(dim * dim)));

    for (uint i_el = 0; i_el < n_iga_elements; ++i_el)
    {
      for (uint i_Pt = 0; i_Pt < n_points_iga_element; ++i_Pt)
      {
        auto& data_scalar_el_Pt = data_scalar[i_el][i_Pt];
        auto& data_vector_el_Pt = data_vector[i_el][i_Pt];
        auto& data_tensor_el_Pt = data_tensor[i_el][i_Pt];

        // Scalar field
        data_scalar_el_Pt[0] = i_Pt % 2;

        // Vector field
        for (uint i = 0; i < dim; ++i)
          data_vector_el_Pt[i] = i_Pt % 2 + i;

        // Tensor field
        for (uint i = 0; i < dim * dim; ++i)
          data_tensor_el_Pt[i] = i_Pt % 2 + i;

      }
    }

    writer.add_point_data(1, "scalar", data_scalar, "scalar field");
    writer.add_point_data(dim, "vector", data_vector, "vector field");
    writer.add_point_data(dim * dim, "tensor", data_tensor, "tensor field");

    string filename = "grid_field" + to_string(dim);
    writer.save(filename);
    writer.save(filename,"appended");
    writer.print_info(out);
}


int main()
{
    test<1>();
    test<2>();
    test<3>();

    return 0;
}

