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
 *  Test for the SphereGridFunction class as a domain
 *
 *  author: martinelli
 *  date: Nov 06, 2015
 *
 */

#include "domain_values.h"


template <int dim>
std::shared_ptr<const Domain<dim,1> >
create_sphere_domain()
{
  auto grid = Grid<dim>::const_create();

  using Sph = grid_functions::SphereGridFunction<dim>;
  return Domain<dim,1>::const_create(Sph::const_create(grid));
}

int main()
{
  out.depth_console(10);

  out.begin_item("SphereGridFunction<1>");
  domain_values<1,1>(*create_sphere_domain<1>(),QUniform<1>::create(3));
  out.end_item();

  out.begin_item("SphereGridFunction<2>");
  domain_values<2,1>(*create_sphere_domain<2>(),QUniform<2>::create(3));
  out.end_item();

  return 0;
}
