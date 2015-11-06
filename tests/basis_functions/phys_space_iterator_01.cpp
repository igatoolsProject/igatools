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

/**
 *  Physical space element evaluation
 *  Function map is the identity function
 *
 *  author: pauletti
 *  date:
 *
 */

#include <igatools/functions/identity_function.h>

#include "phys_space_iterator.h"

template <int dim>
std::shared_ptr<GridFunction<dim,dim>>
                                    create_identity_function(std::shared_ptr<Grid<dim>> &grid)
{
  using Function = grid_functions::LinearGridFunction<dim,dim>;
  typename Function::Value    b;
  typename Function::Gradient A;
  for (int i=0; i<dim; i++)
  {
    A[i][i] = 1.;
    b[i] = 0.0;
  }

  return Function::create(grid,A,b);
}

template<int dim, int sub_dim = dim>
void
identity_map(const int n_knots, const int deg, const string prop=DofProperties::active,
             const bool use_bdry=true)
{
  OUTSTART

  using std::to_string;
  out << "identity_map<" << to_string(dim) << "," << to_string(sub_dim)
      << ">(" << to_string(n_knots) << "," << to_string(deg)
      << "," <<  prop << "," << use_bdry << ")" << std::endl;

  auto grid  = Grid<dim>::create(n_knots);
  auto grid_func = create_identity_function(grid);

  auto phys_basis = create_phys_basis<dim>(grid, grid_func, deg);
  /*
  out.begin_item("Basis");
  phys_basis->print_info(out);
  out.end_item();
  //*/
  const int n_qp = 1;
  elem_values<dim, sub_dim>(phys_basis, n_qp, prop, use_bdry);

  OUTEND
}


template <int dim>
void
dim_test()
{
  out << "dim_test<" << std::to_string(dim) << ">()" << std::endl;
//  identity_map<dim,dim>(3,2);

  for (int n_knots=2; n_knots<4; ++n_knots)
    for (int deg=1; deg<3; ++deg)
    {
      identity_map<dim,dim>(n_knots, deg);
      identity_map<dim, dim-1>(n_knots, deg);
    }
  //*/
}


template<int dim, int sub_dim = dim>
void
identity_map_prop(const int n_knots, const int deg, const bool use_bdry=true)
{
  OUTSTART

  using std::to_string;
  out << "identity_map_prop<" << to_string(dim) << "," << to_string(sub_dim)
      << ">(" << to_string(n_knots) << "," << to_string(deg)
      << "," <<  use_bdry << ")" << std::endl;

  auto grid  = Grid<dim>::create(n_knots);
  auto grid_func = create_identity_function(grid);
  auto space = create_space_prop<dim>(grid, grid_func, deg);
  const int n_qp = 1;
  elem_values<dim, sub_dim>(space, n_qp, DofProp::interior, use_bdry);
  elem_values<dim, sub_dim>(space, n_qp, DofProp::dirichlet, use_bdry);
  elem_values<dim, sub_dim>(space, n_qp, DofProp::neumman, use_bdry);
  OUTEND
}


int main()
{
  dim_test<1>();
  dim_test<2>();
  dim_test<3>();

  identity_map_prop<1,1>(2,1,true);
  identity_map_prop<1,0>(3,1,false);
  identity_map_prop<2,1>(3,2,false);
//*/
  return 0;
}
