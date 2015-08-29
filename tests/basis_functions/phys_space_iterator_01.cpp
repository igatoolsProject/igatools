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


template<int dim, int sub_dim = dim>
void
identity_map(const int n_knots, const int deg, const string prop=DofProperties::active,
             const bool use_bdry=true)
{
  OUTSTART

  auto grid  = CartesianGrid<dim>::create(n_knots);
  auto map = IdentityFunction<dim>::create(grid);
  auto space = create_space<dim>(grid, map, deg);
  const int n_qp = 1;
  elem_values<dim, sub_dim>(space, n_qp, prop, use_bdry);

  OUTEND
}


template <int dim>
void
dim_test()
{
  for (int n_knots=2; n_knots<4; ++n_knots)
    for (int deg=1; deg<3; ++deg)
    {
      identity_map<dim,dim>(n_knots, deg);
      identity_map<dim, dim-1>(n_knots, deg);
    }
}


template<int dim, int sub_dim = dim>
void
identity_map_prop(const int n_knots, const int deg, const bool use_bdry=true)
{
  OUTSTART

  auto grid  = CartesianGrid<dim>::create(n_knots);
  auto map = IdentityFunction<dim>::create(grid);
  auto space = create_space_prop<dim>(grid, map, deg);
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

  return 0;
}
