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
 * Testing the writer, this the add_field()
 * author: martinelli
 * date: Nov 25, 2014
 *
 */

#include "../tests.h"
#include "igatools/io/writer.h"
#include "igatools/functions/grid_function_lib.h"
#include "igatools/functions/function_lib.h"



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


template<int dim>
void
test()
{
  OUTSTART

  const int n_knots = 4;
  auto grid = Grid<dim>::const_create(n_knots);

  auto domain = create_domain_from_grid<dim,0>(grid);

  Writer<dim> writer(domain,2);

//  auto identity_function = IdentityFunction<dim>::create(grid);

  using ScalarFunc = functions::ConstantFunction<dim,0,1,1>;
  using ScalarValue = typename ScalarFunc::Value;
  ScalarValue scalar_value({1.0});
  auto scalar_function = ScalarFunc::const_create(domain,scalar_value);

  using VectorFunc = functions::ConstantFunction<dim,0,dim,1>;
  using VectorValue = typename VectorFunc::Value;
  VectorValue vector_value;
  for (int i = 0 ; i < dim ; ++i)
    vector_value[i] = i;
  auto vector_function = VectorFunc::const_create(domain,vector_value);

#if 0
  using TensorFunc = functions::ConstantFunction<dim,0,dim,2>;
  using TensorValue = typename TensorFunc::Value;
  TensorValue tensor_value;
  for (int i = 0 ; i < dim ; ++i)
    for (int j = 0 ; j < dim ; ++j)
      tensor_value[i][j] = i*dim + j;
  shared_ptr<const Function<dim,0,dim,2>> tensor_function = TensorFunc::create(grid,identity_function,tensor_value);
#endif

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

  string filename = "grid" + to_string(dim);
  writer.save(filename, "ascii");
  writer.save(filename + "_bin", "appended");
  writer.print_info(out);

  out << endl;

  OUTEND
}


int main()
{
  test<1>();
  test<2>();
  test<3>();

  return 0;
}

