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



void eig_dense_matrix(const DenseMatrix &A,
                      SafeSTLVector<Real> &eigenvalues_real,
                      SafeSTLVector<Real> &eigenvalues_imag,
                      DenseMatrix &eigenvectors)
{
    const int n_rows = A.get_num_rows();
    const int n_cols = A.get_num_cols();
    Assert(n_rows == n_cols,ExcDimensionMismatch(n_rows,n_cols));

    const int n = n_rows;

    Assert(eigenvalues_real.size() == n,ExcDimensionMismatch(eigenvalues_real.size(),n));
    Assert(eigenvalues_imag.size() == n,ExcDimensionMismatch(eigenvalues_imag.size(),n));

    Assert(eigenvectors.get_num_rows() == n,ExcDimensionMismatch(eigenvectors.get_num_rows(),n));
    Assert(eigenvectors.get_num_cols() == n,ExcDimensionMismatch(eigenvectors.get_num_cols(),n));

    DenseMatrix eigenvectors_trans(n,n);

    Teuchos::LAPACK<int, double> lapack;
    int info;

    const int workspace_size = 4*n; // thsi should be >= 3*n (higher values gives better performances)
    double workspace[workspace_size];

    double dummy;

    // The Lapack GEEV routine assumes that the matrix entries are sorted column-wise
    // (i.e. using the Fortran way), but in C/C++ the matrix entries are sorted row-wise,
    // so we need to work on the transposte of the input matrix
    DenseMatrix A_trans = boost::numeric::ublas::trans(A);

    /*
    out.begin_item("Matrix A_trans:");
    A_trans.print_info(out);
    out.end_item();
    //*/

    const char jobvl = 'N'; // the left  eigenvector is not computed
    const char jobvr = 'V'; // the right eigenvector is computed
    lapack.GEEV(jobvl,jobvr,n,const_cast<Real *>(&(A_trans.data()[0])),n,
                eigenvalues_real.data(),
                eigenvalues_imag.data(),
                &dummy,1,
                &(eigenvectors_trans.data()[0]),n,
                workspace,
                workspace_size,
                &info);
    if (info > 0)
    {
        AssertThrow(false,ExcMessage("The QR algorithm failed to compute all the "
                                     "eigenvalues, and no eigenvectors have been computed. "
                                     "Elements " + std::to_string(info) + ":" + std::to_string(n-1) +
                                     " of eigenvalues_re and eigenvalues_im contain eigenvalues "
                                     "which have converged."));
    }
    else if (info < 0)
    {
        AssertThrow(false,ExcMessage("The " + std::to_string(std::abs(info)-1) + "-th argument " +
                                     "had an illegal value."));
    }



    std::map<Real,Index> eigval_real_sort;
    for (int i = 0 ; i < n ; ++i)
        eigval_real_sort[eigenvalues_real[i]] = i;


    const auto eigenvalues_imag_tmp = eigenvalues_imag;
    int i = 0;
    for (const auto &eigval_id : eigval_real_sort)
    {
        const int id = eigval_id.second;
//  out << "eigenvalue= " << eigval_id.first <<  "   id=" << eigval_id.second << endl;
        eigenvalues_real[i] = eigval_id.first;
        eigenvalues_imag[i] = eigenvalues_imag_tmp[id];

        // transposing and reordering the eigenvectors
        for (int row = 0 ; row < n ; ++row)
            eigenvectors(row,i) = eigenvectors_trans(id,row);

        ++i;
    }
}


void eig_dense_matrix_symm(const DenseMatrix &A,
                           SafeSTLVector<Real> &eigenvalues,
                           DenseMatrix &eigenvectors)
{
    const int n_rows = A.get_num_rows();
    const int n_cols = A.get_num_cols();
    Assert(n_rows == n_cols,ExcDimensionMismatch(n_rows,n_cols));

    const int n = n_rows;

    Assert(eigenvalues.size() == n,ExcDimensionMismatch(eigenvalues.size(),n));

    Assert(eigenvectors.get_num_rows() == n,ExcDimensionMismatch(eigenvectors.get_num_rows(),n));
    Assert(eigenvectors.get_num_cols() == n,ExcDimensionMismatch(eigenvectors.get_num_cols(),n));

    Teuchos::LAPACK<int, double> lapack;
    int info;

    const int workspace_size = 10*n; // thsi should be >= 3*n (higher values gives better performances)
    double workspace[workspace_size];

    // The Lapack GEEV routine assumes that the matrix entries are sorted column-wise
    // (i.e. using the Fortran way), but in C/C++ the matrix entries are sorted row-wise,
    // so we need to transpose the output matrix
    DenseMatrix eigenvectors_trans = A;
    const char jobz = 'V'; // computes eigenvalues and eigenvectors
    const char uplo = 'L'; // using the lower triangular part of the matrix
    lapack.SYEV(jobz,
                uplo,
                n,
                const_cast<Real *>(&(eigenvectors_trans.data()[0])),
                n,
                eigenvalues.data(),
                workspace,
                workspace_size,
                &info);
    if (info > 0)
    {
        AssertThrow(false,ExcMessage("The algorithm failed to converge. " +
                                     std::to_string(info -1) +
                                     " off-diagonal elements of an intermediate tridiagonal "
                                     "form did not converge to zero."));
    }
    else if (info < 0)
    {
        AssertThrow(false,ExcMessage("The " + std::to_string(std::abs(info)-1) + "-th argument " +
                                     "had an illegal value."));
    }



    std::map<Real,Index> eigval_real_sort;
    for (int i = 0 ; i < n ; ++i)
        eigval_real_sort[eigenvalues[i]] = i;


    int i = 0;
    for (const auto &eigval_id : eigval_real_sort)
    {
        const int id = eigval_id.second;
//  out << "eigenvalue= " << eigval_id.first <<  "   id=" << eigval_id.second << endl;

        eigenvalues[i] = eigval_id.first;

        // transposing and reordering the eigenvectors
        for (int row = 0 ; row < n ; ++row)
            eigenvectors(row,i) = eigenvectors_trans(id,row);
//*/
        ++i;
    }
}

#if 0
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
#endif

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
