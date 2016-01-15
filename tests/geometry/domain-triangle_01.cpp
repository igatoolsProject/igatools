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
 *  Test for a domain built using a TriangleGridFunction
 *
 *  author: martinelli
 *  date: Jan 15, 2016
 */

#if 0
#include "../tests.h"

#include <igatools/geometry/domain_element.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/grid_function_lib.h>
#endif

#include <igatools/io/writer.h>
#include "domain_values.h"

template<int codim>
void test_triangle()
{
  const int dim = 2;

  using std::to_string;
  out.begin_item("test_triangle<codim=" + to_string(codim));

  const int space_dim = dim + codim;
  using F = grid_functions::TriangleGridFunction<space_dim>;

  using Vertex = Points<space_dim>;
  Vertex P0;
  Vertex P1;
  Vertex P2;

  if (space_dim == 2)
  {
    P0[0] = 0.0;
    P0[1] = 0.0;

    P1[0] = 1.0;
    P1[1] = 0.0;

    P2[0] = 1.0;
    P2[1] = 1.0;
  }
  else if (space_dim == 3)
  {
    P0[0] = 0.0;
    P0[1] = 0.0;
    P0[2] = 0.0;

    P1[0] = 1.0;
    P1[1] = 0.0;
    P1[2] = 0.0;

    P2[0] = 1.0;
    P2[1] = 1.0;
    P2[2] = 1.0;
  }
  else
  {
    AssertThrow(false,ExcNotImplemented());
  }

  auto quad = QGauss<dim>::const_create(2);

  SafeSTLArray<SafeSTLVector<Real>,2> knots{{0.0,1.0,2.0},{0.0,0.5}};

  auto grid = Grid<dim>::create(knots);
  auto grid_func = F::create(grid,P0,P1,P2);
  auto domain = Domain<dim,codim>::create(grid_func);


  domain_values<dim,codim>(*domain,quad);

  std::string filename_grid = "grid_" + std::to_string(space_dim) + "D";
  Writer<dim,codim> writer_grid(grid,2);
  writer_grid.save(filename_grid);

  std::string filename_domain = "triangle_" + std::to_string(space_dim) + "D";
  Writer<dim,codim> writer_domain(domain,2);
  writer_domain.save(filename_domain);

  out.end_item();
}


int main()
{
  test_triangle<0>();
  test_triangle<1>();

  return 0;
}

