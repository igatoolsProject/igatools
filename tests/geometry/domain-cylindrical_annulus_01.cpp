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
 *  Test cylindrical annulus mapping values.
 *  author: martinelli
 *  date: Nov 06, 2015
 *
 */


#include "domain_values.h"


std::shared_ptr<const GridFunction<3,3> >
create_cylinder_function()
{
  auto grid = Grid<3>::create({{1.,2.},{0.,numbers::PI/2.},{0.,2.}});

  using Cyl = grid_functions::CylindricalAnnulusGridFunction;
  auto cylinder = Cyl::create(grid);

  return cylinder;
}

int main()
{
  out.depth_console(10);

  auto cyl = create_cylinder_function();
  auto domain = Domain<3>::const_create(cyl);

  auto quad = QGauss<3>::create(1);

  domain_values<3,0>(*domain,quad);

  return 0;
}


