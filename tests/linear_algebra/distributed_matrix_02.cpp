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
 *  Test for definition of non-square Epetra Matrices using sparsity pattern
 *  created by using two different spaces.
 *  this test i going to be an adaptation of:
 *  matrix_definition01.cpp
 *  author: antolin
 *  date: 2013-04-10
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/linear_algebra/epetra_matrix.h>

using namespace EpetraTools;

template<int dim = 1, int range  = 1, int rank=1>
void fill_matrix_and_vector()
{
  using Space = BSplineSpace<dim, range, rank>;
  const int p_r = 3;
  const int p_c = 2;

  out << " Domain dim: " << dim;
  out << " Range dim: " << range <<endl;
  out << " Degree of rows space: " << p_r <<endl;
  out << " Degree of columns space: " << p_c <<endl;

  auto grid = CartesianGrid<dim>::create();
  auto r_space = Space::create(p_r, grid);
  auto c_space = Space::create(p_c, grid);

  const auto n_basis_sp_rows = r_space->get_num_basis();
  const auto n_basis_sp_cols = c_space->get_num_basis();
  out << endl;
  out << "Number of dofs of rows space: " << n_basis_sp_rows << endl;
  out << "Number of dofs of columns space: " << n_basis_sp_cols << endl;
  out << endl;

  Epetra_SerialComm comm;
  auto graph = create_graph(*r_space, DofProperties::active,
                            *c_space, DofProperties::active, comm);


  auto A = create_matrix(*graph);
  A->FillComplete();

  out << "A matrix" << endl;
  A->print_info(out);
  out << endl;
}


int main()
{
  fill_matrix_and_vector();
  return 0;
}

