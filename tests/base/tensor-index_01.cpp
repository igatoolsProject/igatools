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
 *  Test for using the linear and tensor indices.
 *
 *  author: pauletti
 *  date: Feb 18, 2013
 *
 */

#include "../tests.h"

#include <igatools/base/tensor.h>

template <int dim>
void run_test0()
{
  out << "Test for printing a TensorIndex<" << dim << ">" << endl;

  TensorIndex<dim> index;

  out << index << endl;
  out << endl;
}


//Test linear index access
template <int dim, class type>
void run_test1()
{
  const int rank = 2;
  out << "Test accesing  a tensor using tensor indices" << endl;
  out << "Tensor of rank: " << rank << " and dim: " << dim << endl;

  Tensor< dim, rank, type, Tdouble > A;

  TensorIndex<dim> index;

  for (int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j)
      A[i][j] = i*dim + j;
  out << A << endl;


  for (int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j)
    {
      index[0] = i;
      index[1] = j;
      out << A(index) << "  ";
      out << A[i][j]  << endl;
    }
  out << endl;
}



template <int dim, class type>
void run_test2()
{
  const int rank = 2;
  out << "Test the access to a tensor using tensor indices," << endl;
  out << "on a tensor of rank: " << rank << " and dim: " << dim << ".\n";

  Tensor< dim, rank, type, Tdouble > A;

  TensorIndex<rank> index;

  for (int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j)
      A[i][j] = i*dim + j;
  out << A << endl;


  const int size = Tensor< dim, rank, type, Tdouble >::size;

  out << "using linear indices" << endl;
  for (int i=0; i<size; ++i)
    out << A(i) << endl;

  out << "transforming linear to tensor indices" << endl;
  for (int i=0; i<size; ++i)
  {
    auto index = A.flat_to_tensor_index(i);
    out << A(i) << "\t" << A(index) << endl;
  }

  out << endl;
}



template <int dim, class type>
void run_test3()
{
  const int rank = 3;
  out << "Test the access to a tensor using tensor indices," << endl;
  out << "on a tensor of rank: " << rank << " and dim: " << dim << ".\n";

  Tensor< dim, rank, type, Tdouble > A;

  TensorIndex<rank> index;

  for (int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j)
      for (int k=0; k<dim; ++k)
        A[i][j][k] = i*(dim*dim) + j*dim + k;
  out << A << endl;


  const int size = Tensor< dim, rank, type, Tdouble >::size;

  out << "using linear indices" << endl;
  for (int i=0; i<size; ++i)
    out << A(i) << endl;

  out << "transforming linear to tensor indices" << endl;
  for (int i=0; i<size; ++i)
  {
    auto index = A.flat_to_tensor_index(i);
    out << A(i) << "\t" << A(index) << endl;
  }

  out << endl;
}


int main()
{
  out.depth_console(10);

  run_test0<0>();
  run_test0<1>();
  run_test0<2>();
  run_test0<3>();

  run_test3<2, tensor::covariant>();
  run_test2<3, tensor::covariant>();

  return  0;
}

