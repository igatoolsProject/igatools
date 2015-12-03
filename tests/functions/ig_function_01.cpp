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
 *  Test for Function class, as a prototype for an spline function
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/functions/ig_function.h>
#include <igatools/geometry/grid_function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>


#include "function_test.h"

template<int dim, int range>
void test()
{
  auto grid = Grid<dim>::create(3);
  const int deg = 2;
  auto space = SplineSpace<dim,range>::create(deg, grid);
  auto ref_basis = BSpline<dim,range>::create(space);

  auto grid_func = grid_functions::IdentityGridFunction<dim>::create(grid);
  auto domain = Domain<dim>::create(grid_func);

  auto phys_basis = PhysicalSpaceBasis<dim,range>::create(ref_basis,domain);


  Epetra_SerialComm comm;
  auto map = EpetraTools::create_map(*phys_basis, "active", comm);
  auto coeff = EpetraTools::create_vector(*map);
  (*coeff)[0] = 1.;

  auto func = IgFunction<dim,0,range,1>::create(phys_basis, *coeff);

  function_values(*func);

#if 0
  auto func_handler = func->create_cache_handler();

  using Flags = function_element::Flags;
  auto flag = Flags::D0 | Flags::D1 | Flags::D2;
  func_handler->set_element_flags(flag);

  auto elem = func->begin();
  auto end  = func->end();

  func_handler->init_cache(*elem,quad);

  for (; elem != end; ++elem)
  {
    func_handler->fill_element_cache(*elem);

    out.begin_item("FunctionElement:");
    elem->get_index().print_info(out);


    out.begin_item("IgFunction Values:");
    elem->template get_values_from_cache<function_element::_D<0>, dim>(0).print_info(out);
    out.end_item();

    out.begin_item("IgFunction Gradients:");
    elem->template get_values_from_cache<function_element::_D<1>, dim>(0).print_info(out);
    out.end_item();

    out.begin_item("IgFunction Hessians:");
    elem->template get_values_from_cache<function_element::_D<2>, dim>(0).print_info(out);
    out.end_item();

    out.end_item();
  }
#endif
}


int main()
{
  test<2,1>();
//    test<3,3>();

  return 0;
}

