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
 *  Test for Function class, as a prototype for an spline function
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/functions/ig_grid_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>


template<int dim, int range>
void test()
{
  auto quad = QGauss<dim>::const_create(2);
  auto grid = Grid<dim>::const_create(3);
  const int deg = 1;
  auto space = SplineSpace<dim,range>::const_create(deg, grid);
  auto basis = BSpline<dim,range>::const_create(space);

  using Function = IgGridFunction<dim,range>;


  Epetra_SerialComm comm;
  auto map = EpetraTools::create_map(*basis, "active", comm);
  auto coeff = EpetraTools::create_vector(*map);
  (*coeff)[0] = 1.;

  auto func = Function::const_create(basis, *coeff);

  auto func_handler = func->create_cache_handler();

  using Flags = grid_function_element::Flags;
  auto flag = Flags::D0 | Flags::D1 | Flags::D2;
  func_handler->set_element_flags(flag);

  auto elem = func->begin();
  auto end  = func->end();

  func_handler->init_cache(*elem,quad);

  for (; elem != end; ++elem)
  {
    func_handler->fill_element_cache(*elem);
//        elem->get_points().print_info(out);
//        out << endl;
    elem->template get_values_from_cache<grid_function_element::_D<0>, dim>(0).print_info(out);
    out << endl;
    elem->template get_values_from_cache<grid_function_element::_D<1>, dim>(0).print_info(out);
    out << endl;
    elem->template get_values_from_cache<grid_function_element::_D<2>, dim>(0).print_info(out);
    out << endl;
  }

}


int main()
{
  test<2,1>();
//    test<3,3>();

  return 0;
}

