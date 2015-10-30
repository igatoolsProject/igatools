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
 * Testing the writer, this the add_field() using an
 * IgFunction built with a BSplineSpace
 * author: martinelli
 * date: Oct 21, 2015
 *
 */

#include "../tests.h"
#include "igatools/io/writer.h"
#include "igatools/geometry/grid_function_lib.h"
#include "igatools/basis_functions/bspline_space.h"
#include "igatools/functions/ig_function.h"
//#include "igatools/functions/function_lib.h"



template<int dim,int codim>
inline
std::shared_ptr<const Domain<dim,codim> >
create_domain_from_grid(const shared_ptr<const Grid<dim>> &grid)
{
  const int space_dim = dim+codim;
  using F = grid_functions::LinearGridFunction<dim,space_dim>;

  using Grad = typename F::Gradient;
  using Val = typename F::Value;

  Grad A;
  for (int i = 0 ; i < dim ; ++i)
    A[i][i] = Tdouble(1.0);

  Val b;

  auto domain = Domain<dim,codim>::const_create(F::const_create(grid,A,b));
  return domain;
}


template<int dim,int codim,int range>
std::shared_ptr<const PhysicalSpace<dim,range,1,codim> >
create_phys_space(const shared_ptr<const Domain<dim,codim>> &domain)
{
  const int deg = 2;

  const auto grid = domain->get_grid_function()->get_grid();

  auto bsp_space = BSplineSpace<dim,range,1>::const_create(deg,grid);

  auto phys_space = PhysicalSpace<dim,range,1,codim>::const_create(bsp_space,domain);

  return phys_space;
}


template<int dim,int codim,int range>
std::shared_ptr<const Function<dim,codim,range,1> >
create_ig_function(const shared_ptr<const Domain<dim,codim>> &domain)
{
  auto phys_space = create_phys_space<dim,codim,range>(domain);

  int n_basis = phys_space->get_num_basis();

  IgCoefficients coeffs;
  for (int dof = 0 ; dof < n_basis ; ++dof)
    coeffs[dof] = 1.0;

  auto ig_func = IgFunction<dim,codim,range,1>::const_create(phys_space,coeffs);

  return ig_func;
}

template<int dim>
void
test()
{
  OUTSTART

  out.begin_item("test<" + std::to_string(dim) + ">");

  const int n_knots = 4;
  auto grid = Grid<dim>::const_create(n_knots);

  auto domain = create_domain_from_grid<dim,0>(grid);

  Writer<dim> writer(domain,2);

  auto scalar_function = create_ig_function<dim,0,1>(domain);

  auto vector_function = create_ig_function<dim,0,dim>(domain);


  SafeSTLVector<Real> cell_data(grid->get_num_all_elems());
  int elem_id = 0;
  int n=1;
  for (const auto &elem : *grid)
  {
    cell_data[elem_id++] = n;
    n *= -1;
  }
  writer.add_element_data(cell_data, "chess board");

  writer.add_field(*scalar_function,"scalar_function");
  writer.add_field(*vector_function,"vector_function");
//    writer.add_field(tensor_function,"tensor_function");

  string filename = "grid_dim" + to_string(dim);
  writer.save(filename);
  writer.save(filename,"appended");
  writer.print_info(out);

  out.end_item();

  OUTEND
}


int main()
{
  test<1>();
  test<2>();
  test<3>();

  return 0;
}

