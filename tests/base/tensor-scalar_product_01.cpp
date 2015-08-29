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
 * Test for scalar product of rank 1 tensors
 * Author: pauletti 2012/10/26
 *
 */
#include "../tests.h"

#include <igatools/base/tensor.h>


//template< int dim >
//void do_scalar_product0()
//{
//  const int rank = 0;
//  out << "Test for Tensor of rank: " << rank << " and dim: " << dim << endl;
//  Tzero< tensor::contravariant, Tdouble > t1, t2;
//  t1(0) = 2;
//  t2(0) = 3;
//  out << "The scalar product of:" << std::endl;
//  out << t1 << std::endl;
//  out << "with:" << std::endl;
//  out << t2 << std::endl;
//  out << "is:" << std::endl;
//  out << scalar_product (t1, t2) << std::endl;
//  out << endl;
//}


template< int dim >
void do_scalar_product1()
{
  const int rank = 1;
  out << "Test for Tensor of rank: " << rank << " and dim: " << dim << endl;
  Tensor < dim, rank, tensor::contravariant, Tdouble > t1, t2;
  for (int i = 0 ; i < dim ; i++)
  {
    t1[i]=i;
    t2[i]=dim+i;
  }

  out << "The scalar product of:" << std::endl;
  out << t1 << std::endl;
  out << "with:" << std::endl;
  out << t2 << std::endl;
  out << "is:" << std::endl;
  out << scalar_product(t1, t2) << std::endl;
  out << endl;
}


template< int dim >
void do_scalar_product2()
{
  const int rank = 2;
  out << "Test for Tensor of rank: " << rank << " and dim: " << dim << endl;
  Tensor < dim, rank, tensor::contravariant, Tdouble > t1, t2;
  for (int i = 0 ; i < dim ; ++i)
    for (int j = 0 ; j < dim ; ++j)
    {
      t1[i][j]=i*dim+j;
      t2[i][j]=pow(dim,rank) + t1[i][j];
    }

  out << "The scalar product of:" << std::endl;
  out << t1 << std::endl;
  out << "with:" << std::endl;
  out << t2 << std::endl;
  out << "is:" << std::endl;
  out << scalar_product(t1, t2) << std::endl;

  out << endl;
}


//template< int dim >
//void do_scalar_product01()
//{
//  const int rank = 1;
//  out << "Test for Tensor of rank: " << rank << " and dim: " << dim << endl;
//  Tzero<tensor::covariant,
//  Tensor < dim, rank, tensor::contravariant, Tdouble >> t1, t2;
//      for ( int i = 0 ; i < dim ; i++ )
//      {
//         t1(0)[i]=i;
//         t2(0)[i]=dim+i;
//      }
//
//      out << "The scalar product of:" << std::endl;
//      out << t1 << std::endl;
//      out << "with:" << std::endl;
//      out << t2 << std::endl;
//      out << "is:" << std::endl;
//      out << scalar_product (t1, t2) << std::endl;
//      out << endl;
//}


template< int dim >
void do_scalar_product11()
{
  const int rank = 1;
  out << "Test for Tensor of rank: " << 2*rank << " and dim: " << dim << endl;
  Tensor < dim, rank, tensor::covariant,
         Tensor < dim, 1, tensor::contravariant, Tdouble >> t1, t2;

  for (int i = 0 ; i < dim ; i++)
    for (int j = 0 ; j < dim ; j++)
    {
      t1[i][j]=i*dim+j;
      t2[i][j]=pow(dim,rank+1) + t1[i][j];
    }

  out << "The scalar product of:" << std::endl;
  out << t1 << std::endl;
  out << "with:" << std::endl;
  out << t2 << std::endl;
  out << "is:" << std::endl;
  out << scalar_product(t1, t2) << std::endl;
  out << endl;
}


int main(int argc, char *argv[])
{

  out.depth_console(10);
//  do_scalar_product0<1>();
//  do_scalar_product0<3>();

  do_scalar_product1<1>();
  do_scalar_product1<2>();
  do_scalar_product1<3>();

  do_scalar_product2<1>();
  do_scalar_product2<2>();
  do_scalar_product2<3>();

//  do_scalar_product01<1>();
//  do_scalar_product01<2>();
//  do_scalar_product01<3>();

  do_scalar_product11<1>();
  do_scalar_product11<2>();
  do_scalar_product11<3>();

}

