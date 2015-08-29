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
 *  Test for the co transpose of a tensor.
 *
 *  author: pauletti
 *  date: apr 18, 2013
 *
 */

#include "../tests.h"

#include <igatools/base/tensor.h>



template <int dim, int range, int rank>
void run_test()
{

  out << "dim: " << dim  << " range: " << range << endl;
  Derivatives<dim, range, rank, 1> DF;
  for (int i = 0; i < dim; ++i)
    DF[i][i] = double(i+1) ;

  Real det;
  auto DF_inv = inverse(DF, det);

  auto DF_inv_t = co_tensor(transpose(DF_inv));

  Points<dim> n_hat;
  for (int i = 0; i < dim; ++i)
    n_hat[i] = double(i+1) ;
  auto n = action(DF,n_hat);

  out << "Action of: " << DF << "on:" << n_hat;
  out << "is:" << n << endl;

  auto n1 = action(DF_inv_t, n_hat);

  out << "Action of: " << DF_inv_t << "on:" << n_hat;
  out << "is:" << n1 << endl;

}


int main()
{
  out.depth_console(10);

  run_test<1,1,1>();
  run_test<2,2,1>();
  run_test<3,3,1>();

  run_test<1,2,1>();
  run_test<1,3,1>();
  run_test<2,3,1>();

  return  0;
}

