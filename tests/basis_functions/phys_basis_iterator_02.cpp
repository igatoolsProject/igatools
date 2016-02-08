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
 *  Test for the evaluation of physical basis functions
 *  with the ball function as a map.
 *
 *  author: pauletti
 *  date: 2014/11/10
 *
 */

#include <igatools/functions/identity_function.h>
#include <igatools/functions/function_lib.h>

#include "phys_basis_iterator.h"



template<int dim, int sub_dim = dim>
void
ball_map(const int n_knots, const int deg, const string prop=DofProperties::active,
         const bool use_bdry=true)
{
  OUTSTART

  BBox<dim> box;
  box[0] = {0.5, 1};
  for (int i=1; i<dim; ++i)
    box[i] = {0., M_PI/4.};

  auto grid  = Grid<dim>::create(box, n_knots);
  auto map = grid_functions::BallGridFunction<dim>::create(grid);
  auto phys_basis = create_phys_basis<dim>(grid, map, deg);
  const int n_qp = 2;
  elem_values<dim, sub_dim>(phys_basis, n_qp, prop, use_bdry);

  OUTEND
}

template<int dim, int sub_dim = dim>
void
ball_map_prop(const int n_knots, const int deg, const bool use_bdry=true)
{
  OUTSTART

  BBox<dim> box;
  box[0] = {0.5, 1};
  for (int i=1; i<dim; ++i)
    box[i] = {0., M_PI/4.};

  auto grid  = Grid<dim>::const_create(box, n_knots);

  auto map = grid_functions::BallGridFunction<dim>::const_create(grid);
  auto phys_basis = create_phys_basis<dim>(grid, map, deg);
  const int n_qp = 1;
  elem_values<dim, sub_dim>(phys_basis, n_qp, DofProp::interior, use_bdry);
  elem_values<dim, sub_dim>(phys_basis, n_qp, DofProp::dirichlet, use_bdry);
  elem_values<dim, sub_dim>(phys_basis, n_qp, DofProp::neumman, use_bdry);
  OUTEND
}






int main()
{

  //  ball_map_prop<1,1>(2,1,true);


  ball_map<1,1>(2,1);
  ball_map<2,1>(3,2);
  ball_map<3,3>(2,1);

  //ball_map_prop<2,1>(3,1);
  //ball_map_prop<2,1>(3,2, false);
  return 0;
}
