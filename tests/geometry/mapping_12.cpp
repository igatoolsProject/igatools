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
void linear_grid_func()
{
  out.begin_item("linear_grid_func<dim=" + std::to_string(dim) +",codim=" + std::to_string(codim) +">");

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
  int elem_id = 0;
  for (; elem != end; ++elem, ++elem_id)
  {
    out.begin_item("Element " + std::to_string(elem_id));

    out << "Element ID: " << elem->get_index() << std::endl;

    domain_handler->fill_element_cache(elem);

    out.begin_item("Inverse Jacobian:");
    elem->template get_values_from_cache<domain_element::_InvJacobian,dim>(0).print_info(out);
    out.end_item();

    out.begin_item("Inverse Hessian:");
    elem->template get_values_from_cache<domain_element::_InvHessian,dim>(0).print_info(out);
    out.end_item();

    out.end_item();
  }

  out.end_item();
}


int main()
{
  linear_grid_func<1,0>();
  linear_grid_func<2,0>();
  linear_grid_func<3,0>();

  linear_grid_func<1,1>();
  linear_grid_func<2,1>();

  linear_grid_func<1,2>();

  return 0;
}

