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
 *  Test for the DenseMatrix eigenvalues
 *
 *  author: pauletti
 *  date: 2014-10-24
 *
 */

#include "../tests.h"
#include <igatools/linear_algebra/dense_matrix.h>

#include <igatools/base/tensor.h>

#include "Teuchos_LAPACK.hpp"


DenseMatrix create_matrix(const int n)
{

  DenseMatrix matrix_eigenvalues(n,n) ;
  matrix_eigenvalues = 0.0;
  for (int i = 0 ; i < n ; ++i)
    matrix_eigenvalues(i,i) = Real(i);

  /*
  out.begin_item("Matrix eigenvalues:");
  matrix_eigenvalues.print_info(out);
  out.end_item();
  //*/

  DenseMatrix matrix_eigenvectors(n,n);
  matrix_eigenvectors = 0.0;
  for (int i = 0 ; i < n ; ++i)
    for (int j = i ; j < n ; ++j)
      matrix_eigenvectors(i,j) = 1.0;
  /*
      out.begin_item("Matrix eigenvectors:");
      matrix_eigenvectors.print_info(out);
      out.end_item();
  //*/

  Real det = 0.0;
  DenseMatrix inverse_matrix_eigenvector = matrix_eigenvectors.inverse(det);

  /*
  out.begin_item("Inverse of matrix eigenvectors:");
  inverse_matrix_eigenvector.print_info(out);
  out.end_item();
  //*/

  using boost::numeric::ublas::prod;
  const DenseMatrix tmp = prod(matrix_eigenvectors,matrix_eigenvalues);
  DenseMatrix A = prod(tmp,inverse_matrix_eigenvector);

  return A ;
}




void do_test_nonsymmetric_matrix(const int n)
{
  out.begin_item("do_test_nonsymmetric_matrix(" + std::to_string(n) + ")");

  DenseMatrix A = create_matrix(n);

  out.begin_item("Matrix A:");
  A.print_info(out);
  out.end_item();

  SafeSTLVector<Real> eigenvalues_real(n);
  SafeSTLVector<Real> eigenvalues_imag(n);
  DenseMatrix eigenvectors(n,n);
  eig_dense_matrix(A,eigenvalues_real,eigenvalues_imag,eigenvectors);

  out.begin_item("eigenvalues real:");
  eigenvalues_real.print_info(out);
  out.end_item();

  out.begin_item("eigenvalues imag:");
  eigenvalues_imag.print_info(out);
  out.end_item();

  out.begin_item("eigenvectors:");
  eigenvectors.print_info(out);
  out.end_item();

  out.end_item();
}


void do_test_symmetric_matrix(const int n)
{
  out.begin_item("do_test_symmetric_matrix(" + std::to_string(n) + ")");

  DenseMatrix tmp = create_matrix(n);
  DenseMatrix A = boost::numeric::ublas::prod(tmp,boost::numeric::ublas::trans(tmp)) ;

  out.begin_item("Matrix A:");
  A.print_info(out);
  out.end_item();

  SafeSTLVector<Real> eigenvalues(n);
  DenseMatrix eigenvectors(n,n);
  eig_dense_matrix_symm(A,eigenvalues,eigenvectors);

  out.begin_item("eigenvalues real:");
  eigenvalues.print_info(out);
  out.end_item();

  /*
  out.begin_item("eigenvalues imag:");
  eigenvalues_imag.print_info(out);
  out.end_item();
  //*/

  out.begin_item("eigenvectors:");
  eigenvectors.print_info(out);
  out.end_item();

  out.end_item();
}


int main()
{
  /*
    eigen_values<2>();
    eigen_values<3>();

    eigen_values2();
  //*/

  const int n_max = 5;
  for (int n = 2 ; n <= n_max ; ++n)
  {
    do_test_nonsymmetric_matrix(n);
    do_test_symmetric_matrix(n);
  }

  return 0;
}
