//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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
 *  Test for building a matrix on a basis of an const igfunction
 *
 *  author: pauletti
 *  date: 2015-03-17
 *
 */

//TODO (pauletti, Apr 11, 2015): rename this test to something more meaningful
#include "../tests.h"

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/functions/ig_grid_function.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/linear_algebra/epetra_matrix.h>
#include <igatools/linear_algebra/epetra_vector.h>


template<int dim, int sp_dim>
void using_const_basis(shared_ptr<const IgGridFunction<dim,sp_dim>> fun)
{
  OUTSTART

  auto basis = fun->get_basis();
  Epetra_SerialComm comm;
  auto matrix = EpetraTools::create_matrix(*basis,DofProperties::active,comm);

  OUTEND
}

template<int dim, int sp_dim>
void using_const_function(shared_ptr<const GridFunction<dim,sp_dim>> fun)
{
  OUTSTART
  using Func =  grid_functions::ConstantGridFunction<dim, sp_dim>;
  typename Func::Value val {0.};
  auto grid = fun->get_grid();
  auto zero = Func::const_create(grid,val);
  OUTEND
}

int main()
{
  const int dim = 2;
  using Space = SplineSpace<dim>;
  using Basis = BSpline<dim>;

  auto grid = Grid<dim>::const_create(5);
  auto space = Space::const_create(1, grid);
  auto basis = Basis::const_create(space);

  const auto n_basis = basis->get_num_basis();

  IgCoefficients coeffs;
  for (int dof = 0 ; dof < n_basis ; ++dof)
    coeffs[dof] = 1.0;


  using IgGridFunc = IgGridFunction<dim,1>;
  auto fun = dynamic_pointer_cast<const IgGridFunc>(IgGridFunc::const_create(basis,coeffs));

  using_const_basis<2,1>(fun);
  using_const_function<2,1>(fun);

  return 0;
}
