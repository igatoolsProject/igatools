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
 *  Test for physical space using a BSpline as reference space.
 *
 *  author: martinelli
 *  date: Sep 24, 2015
 *
 */

#include "../tests.h"

//#include <igatools/base/quadrature_lib.h>
//#include <igatools/functions/function_lib.h>
//#include <igatools/functions/identity_function.h>

//#include <igatools/basis_functions/physical_basis_element_handler.h>
//#include <igatools/basis_functions/bspline_element.h>
//#include <igatools/basis_functions/physical_basis_element.h>
//#include <igatools/geometry/push_forward_element.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/basis_functions/physical_basis.h>

#if 0
template <int dim, int range=1, int rank=1, int codim = 0>
void cache_init(const ValueFlags flag,
                const int n_knots = 5, const int deg=1)
{
  OUTSTART

  using BspSpace = BSpline<dim, range, rank>;
  using Basis    = PhysicalBasis<dim,range,rank,codim, Transformation::h_grad>;
  auto grid      = Grid<dim>::const_create(n_knots);
  auto ref_space = BspSpace::const_create(deg, grid);

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

  auto quad = QGauss<dim>(2);
  auto linear_func = Function::const_create(grid,A, b);
  auto phys_domain = Domain<dim,codim>::const_create(linear_func);
  auto space = Basis::const_create(ref_space, phys_domain);


  auto elem_handler = space->create_cache_handler();
  elem_handler->reset(flag, quad);
  elem_handler->print_info(out);

  OUTEND
}



template <int dim, int range=1, int rank=1, int codim = 0>
void cache_init_elem(const ValueFlags flag,
                     const int n_knots = 5, const int deg=1)
{
//    const int k = dim;
  OUTSTART

  using BspSpace = BSpline<dim, range, rank>;
  using Basis    = PhysicalBasis<dim,range,rank,codim, Transformation::h_grad>;

  auto grid  = Grid<dim>::const_create(n_knots);
  auto ref_space = BspSpace::const_create(deg, grid);

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

  auto quad = QGauss<dim>(2);
  auto linear_func = Function::const_create(grid,A, b);
  auto phys_domain = Domain<dim,codim>::const_create(linear_func);
  auto space = Basis::const_create(ref_space, phys_domain);

  auto elem_handler = space->create_cache_handler();
  elem_handler->reset(flag, quad);

  auto elem = space->begin();
  elem_handler->init_element_cache(elem);
  elem->print_cache_info(out);

  OUTEND
}


template <int dim, int range=1, int rank=1, int codim = 0>
void cache_fill_elem(const ValueFlags flag,
                     const int n_knots = 5, const int deg=1)
{
  OUTSTART

//   const int k = dim;
  using BspSpace = BSpline<dim, range, rank>;
  using Basis    = PhysicalBasis<dim,range,rank,codim, Transformation::h_grad>;

  auto grid  = Grid<dim>::const_create(n_knots);
  auto ref_space = BspSpace::const_create(deg, grid);

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

  auto quad = QGauss<dim>(2);
  auto linear_func = Function::const_create(grid,A, b);
  auto phys_domain = Domain<dim,codim>::const_create(linear_func);
  auto space = Basis::const_create(ref_space, phys_domain);

  auto elem_handler = space->create_cache_handler();
  elem_handler->reset(flag, quad);

  auto elem = space->begin();
  auto end = space->end();
  elem_handler->init_element_cache(elem);
  for (; elem != end; ++elem)
  {
    elem_handler->fill_element_cache(elem);
    elem->print_cache_info(out);
  }

  OUTEND
}



template <int dim, int range=1, int rank=1, int codim = 0>
void cache_get_elem_values(const ValueFlags flag,
                           const int n_knots = 5, const int deg=1)
{
  OUTSTART
  const int k = dim;
  using BspSpace = BSpline<dim, range, rank>;
  using Basis    = PhysicalBasis<dim,range,rank,codim, Transformation::h_grad>;

  auto grid  = Grid<dim>::const_create(n_knots);
  auto ref_space = BspSpace::const_create(deg, grid);

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

  auto quad = QGauss<dim>(2);
  auto linear_func = Function::const_create(grid,A, b);
  auto phys_domain = Domain<dim,codim>::const_create(linear_func);
  auto space = Basis::const_create(ref_space, phys_domain);

  auto elem_handler = space->create_cache_handler();
  elem_handler->reset(flag, quad);

  auto elem = space->begin();
  auto end = space->end();
  elem_handler->init_element_cache(elem);
  for (; elem != end; ++elem)
  {
    elem_handler->fill_element_cache(elem);
    elem->template get_basis_data<_Value, k>(0,DofProperties::active).print_info(out);
  }

  OUTEND
}
#endif



template <int dim, int range=1, int rank=1, int codim = 0>
std::shared_ptr<const PhysicalBasis<dim,range,rank,codim>>
                                                        create_phys_space()
{
  OUTSTART
  auto grid = Grid<dim>::const_create();
  const int deg = 2;
  auto ref_space = BSpline<dim,range,rank>::const_create(
                     SplineSpace<dim,range,rank>::const_create(deg,grid));

  using GridFunc = grid_functions::BallGridFunction<dim>;
  auto grid_func = GridFunc::const_create(grid);

  using Domain = Domain<dim,codim>;
  auto domain = Domain::const_create(grid_func);

  using PhysSpace = PhysicalBasis<dim,range,rank,codim>;
  auto phys_space = PhysSpace::const_create(ref_space,domain,Transformation::h_grad);


  using std::to_string;
  out.begin_item("PhysicalBasis<"
                 + to_string(dim) + ","
                 + to_string(range) + ","
                 + to_string(rank) + ","
                 + to_string(codim) + ",Transformation::h_grad>");
  phys_space->print_info(out);
  out.end_item();

  OUTEND

  return phys_space;
}

template <int dim, int range=1, int rank=1, int codim = 0>
void
test_phys_space_accessor()
{
  auto phys_space = create_phys_space<dim,range,rank,codim>();

}


int main()
{
  out.depth_console(10);

  test_phys_space_accessor<1>();
//  create_space<2>();
//  create_space<3>();

  return  0;
}
