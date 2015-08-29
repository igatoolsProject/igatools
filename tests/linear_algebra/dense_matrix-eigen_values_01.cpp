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
 *  Test for the DenseMatrix eigenvalues
 *
 *  author: pauletti
 *  date: 2014-10-24
 *
 */

#include "../tests.h"
#include <igatools/linear_algebra/dense_matrix.h>

#include <igatools/base/tensor.h>

template <int dim>
void eigen_values()
{
  OUTSTART

  DenseMatrix A(dim, dim);
  A.clear();
  for (int i=0; i<dim; ++i)
    A(i,i) = i+1;

  A.print_info(out);
  out << endl << "Eigen Values:" << endl;
  A.eigen_values().print_info(out);
  out << endl;

  OUTEND
}


void eigen_values2()
{
  OUTSTART

  const int dim=2;
  Tensor<dim, 1, tensor::covariant, Tensor<dim, 1, tensor::contravariant, Tdouble> > B;
  B[0][0] = 1;
  B[0][1] = 2;
  B[1][0] = 3;
  B[1][1] = 4;
  const auto A = unroll_to_matrix(B);

  out << endl << "Eigen Values:" << endl;
  A.eigen_values().print_info(out);
  out << endl;

  OUTEND
}

int main()
{
  eigen_values<2>();
  eigen_values<3>();

  eigen_values2();

  return 0;
}
