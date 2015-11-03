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
 *  Test for physical space element values
 *
 *  author: pauletti
 *  date: 2014-11-08
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid_function_lib.h>

#include <igatools/basis_functions/phys_space_element_handler.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/physical_space_element.h>



template<int dim, int codim>
using MapFunc = Function<dim, 0, dim+codim>;


template<int dim, int codim=0>
auto
create_function(shared_ptr<const Grid<dim>> grid)
{

  using Function = grid_functions::LinearGridFunction<dim,dim+codim>;
  typename Function::Value    b;
  typename Function::Gradient A;
  for (int i=0; i<dim+codim; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }

  return Function::const_create(grid,A, b);
}


template <int dim, int order = 0, int range=1, int rank=1, int codim = 0>
void elem_values(const int n_knots = 5, const int deg=1)
{
  OUTSTART
  const int k = dim;
  using BspSpace = BSpline<dim, range, rank>;
  using RefSpace = ReferenceSpace<dim, range,rank>;
  using Space = PhysicalSpace<dim,range,rank,codim>;
//    using ElementHandler = typename Space::ElementHandler;

  auto grid  = Grid<dim>::const_create(n_knots);

  auto ref_space = BspSpace::const_create(SplineSpace<dim,range,rank>::const_create(deg,grid));
  auto map_func = create_function(grid);

  auto space = Space::const_create(
                 ref_space,
                 Domain<dim,codim>::const_create(map_func), Transformation::h_grad);

  using space_element::Flags;

  auto flag = Flags::none;
  switch (order)
  {
    case 0:
      flag |= Flags::value;
      break;
    case 1:
      flag |= Flags::gradient;
      break;
    case 2:
      flag |= Flags::hessian;
      break;
  }

  auto quad = QGauss<k>::create(2);

  auto elem_filler = space->create_cache_handler();
  elem_filler->template set_flags<dim>(flag);

  auto elem = space->begin();
  auto end = space->end();
  using space_element::_Value;
  elem_filler->init_element_cache(elem,quad);
  for (; elem != end; ++elem)
  {
    elem_filler->fill_element_cache(elem);
    elem->template get_basis_data<_Value, k>(0,DofProperties::active).print_info(out);
  }

  OUTEND
}



int main()
{
  out.depth_console(10);


  elem_values<1>();

  return  0;
}
