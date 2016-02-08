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
 *  Test for physical basis
 *
 *  author: pauletti
 *  date: Oct 08, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/grid_function_lib.h>

#include <igatools/basis_functions/physical_basis_handler.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/physical_basis_element.h>
//#include <igatools/geometry/push_forward_element.h>


template <int dim, int range=1, int rank=1, int codim = 0>
void cache_init(const basis_element::Flags flag,
                const int n_knots = 5, const int deg=1)
{
  OUTSTART

  using BspBasis = BSpline<dim, range, rank>;
  using Basis    = PhysicalBasis<dim,range,rank,codim>;
  auto grid      = Grid<dim>::const_create(n_knots);
  auto ref_basis = BspBasis::const_create(SplineSpace<dim,range,rank>::const_create(deg,grid));

  using Function = grid_functions::LinearGridFunction<dim,dim+codim>;
  typename Function::Value    b;
  typename Function::Gradient A;
  for (int i=0; i<Basis::space_dim; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }

//  auto quad = QGauss<dim>::create(2);
  auto map_func = Function::const_create(grid,A, b);
  auto basis = Basis::const_create(
                 ref_basis,
                 Domain<dim,codim>::const_create(map_func),Transformation::h_grad);


  auto elem_handler = basis->create_cache_handler();
  elem_handler->template set_flags<dim>(flag);
  elem_handler->print_info(out);

  OUTEND
}



template <int dim, int range=1, int rank=1, int codim = 0>
void cache_init_elem(const basis_element::Flags flag,
                     const int n_knots = 5, const int deg=1)
{
//    const int k = dim;
  OUTSTART

  using BspBasis = BSpline<dim, range, rank>;
  using Basis    = PhysicalBasis<dim,range,rank,codim>;

  auto grid  = Grid<dim>::const_create(n_knots);
  auto ref_basis = BspBasis::const_create(SplineSpace<dim,range,rank>::const_create(deg,grid));

  using Function = grid_functions::LinearGridFunction<dim,dim+codim>;
  typename Function::Value    b;
  typename Function::Gradient A;
  for (int i=0; i<Basis::space_dim; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }

  auto quad = QGauss<dim>::create(2);
  auto map_func = Function::const_create(grid, A, b);
  auto basis = Basis::const_create(
                 ref_basis,
                 Domain<dim,codim>::const_create(map_func), Transformation::h_grad);

  auto elem_handler = basis->create_cache_handler();
  elem_handler->template set_flags<dim>(flag);

  auto elem = basis->begin();
  elem_handler->init_element_cache(elem,quad);
  elem->print_cache_info(out);

  OUTEND
}


template <int dim, int range=1, int rank=1, int codim = 0>
void cache_fill_elem(const basis_element::Flags flag,
                     const int n_knots = 5, const int deg=1)
{
  OUTSTART

//   const int k = dim;
  using BspBasis = BSpline<dim, range, rank>;
  using Basis    = PhysicalBasis<dim,range,rank,codim>;

  auto grid  = Grid<dim>::const_create(n_knots);
  auto ref_basis = BspBasis::const_create(SplineSpace<dim,range,rank>::const_create(deg,grid));

  using Function = grid_functions::LinearGridFunction<dim,dim+codim>;
  typename Function::Value    b;
  typename Function::Gradient A;
  for (int i=0; i<Basis::space_dim; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }

  auto quad = QGauss<dim>::create(2);
  auto map_func = Function::const_create(grid,A, b);
  auto basis = Basis::const_create(
                 ref_basis,
                 Domain<dim,codim>::const_create(map_func), Transformation::h_grad);

  auto elem_handler = basis->create_cache_handler();
  elem_handler->template set_flags<dim>(flag);

  auto elem = basis->begin();
  auto end = basis->end();
  elem_handler->init_element_cache(elem,quad);
  for (; elem != end; ++elem)
  {
    elem_handler->fill_element_cache(elem);
    elem->print_cache_info(out);
  }

  OUTEND
}



template <int dim, int range=1, int rank=1, int codim = 0>
void cache_get_elem_values(const basis_element::Flags flag,
                           const int n_knots = 5, const int deg=1)
{
  OUTSTART
  const int k = dim;
  using BspBasis = BSpline<dim, range, rank>;
  using Basis    = PhysicalBasis<dim,range,rank,codim>;

  auto grid  = Grid<dim>::const_create(n_knots);
  auto ref_basis = BspBasis::const_create(SplineSpace<dim,range,rank>::const_create(deg,grid));

  using Function = grid_functions::LinearGridFunction<dim,dim+codim>;
  typename Function::Value    b;
  typename Function::Gradient A;
  for (int i=0; i<Basis::space_dim; i++)
  {
    for (int j=0; j<dim; j++)
      if (j == i)
        A[j][j] = 2.;
    b[i] = i;
  }

  auto quad = QGauss<dim>::create(2);
  auto map_func = Function::const_create(grid,  A, b);
  auto basis = Basis::const_create(
                 ref_basis,
                 Domain<dim,codim>::const_create(map_func), Transformation::h_grad);

  auto elem_handler = basis->create_cache_handler();
  elem_handler->template set_flags<dim>(flag);

  using basis_element::_Value;

  auto elem = basis->begin();
  auto end = basis->end();
  elem_handler->init_element_cache(elem,quad);
  for (; elem != end; ++elem)
  {
    elem_handler->fill_element_cache(elem);
    elem->template get_basis_data<_Value, k>(0,DofProperties::active).print_info(out);
  }

  OUTEND
}



int main()
{
  out.depth_console(10);

  cache_init<1>(basis_element::Flags::value);
  cache_init_elem<1>(basis_element::Flags::value);
  cache_fill_elem<1>(basis_element::Flags::value);
  cache_get_elem_values<1>(basis_element::Flags::value);

  return  0;
}
