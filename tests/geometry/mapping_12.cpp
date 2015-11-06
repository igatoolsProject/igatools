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
 *  Test for linear mapping class
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/geometry/grid_function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>

template<int dim, int codim>
void test()
{
  const int space_dim = dim + codim;
  using Function = grid_functions::LinearGridFunction<dim,space_dim>;

  typename Function::Value    b;
  typename Function::Gradient A;
  for (int i=0; i<space_dim; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }

  auto grid = Grid<dim>::const_create(3);
  auto F = Function::const_create(grid, A, b);

  using Flags = domain_element::Flags;
  auto flag = Flags::inv_hessian|Flags::inv_jacobian;

  auto domain = Domain<dim,codim>::const_create(F);

  auto domain_handler = domain->create_cache_handler();

  domain_handler->set_element_flags(flag);

  auto elem = domain->begin();
  auto end  = domain->end();


  auto quad = QGauss<dim>::create(2);
  domain_handler->init_cache(elem,quad);
  for (; elem != end; ++elem)
  {
    domain_handler->fill_element_cache(elem);
    out << "Inverse Jacobian:" << std::endl;
    elem->template get_values_from_cache<domain_element::_InvJacobian,dim>(0).print_info(out);
    out << endl;
    out << "Inverse Hessian:" << std::endl;
    elem->template get_values_from_cache<domain_element::_InvHessian,dim>(0).print_info(out);
    out << endl;
  }

}


int main()
{
  test<2,0>();
//    test<3,3>();

  return 0;
}

