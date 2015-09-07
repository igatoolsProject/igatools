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

#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/base/exceptions.h>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "Teuchos_LAPACK.hpp"

#include <algorithm>

IGA_NAMESPACE_OPEN
#if 0
DenseMatrix::
DenseMatrix(const Index n_rows, const Index n_cols, const bool init_to_zero)
  :
  BoostMatrix(n_rows,n_cols)
{
  if (init_to_zero)
    BoostMatrix::clear();
}

DenseMatrix::
DenseMatrix(const Index n,const bool init_to_zero)
  :
  DenseMatrix(n,n,init_to_zero)
{}
#endif

DenseMatrix &
DenseMatrix::operator=(const Real value)
{
  using zero_matrix = boost::numeric::ublas::zero_matrix<Real>;
  Assert(value==0, ExcNonZero());
  *this = zero_matrix(this->size1(), this->size2());
  return *this;
}



auto
DenseMatrix::
get_row(const int row_id) const -> MatrixRowType
{
  return MatrixRowType(*this,row_id);
}



int
DenseMatrix::
size1() const
{
  return int(boost::numeric::ublas::matrix<Real>::size1());
}



int
DenseMatrix::
size2() const
{
  return int(boost::numeric::ublas::matrix<Real>::size2());
}


SafeSTLVector<Real> DenseMatrix::eigen_values() const
{
  Assert(this->size1()==this->size2(), ExcMessage("Should be square"));

  Teuchos::LAPACK<int, double> lapack;
  const int n = this->size1();
  BoostMatrix A(*this);
  SafeSTLVector<Real> e_values_re(n);
  SafeSTLVector<Real> e_values_im(n);
  double e_v[6];
  int info;
  double ws[3*n];
//    lapack.SYEV('V','U', n, &(A.data()[0]), n, &(e_values[0]),
//                ws, 3*n-1, &info);
  lapack.GEEV('N','N', n, &(A.data()[0]), n,
              &(e_values_re[0]), &(e_values_im[0]),
              e_v, n, e_v, n,
              ws, 3*n, &info);

  Assert(info == 0, ExcMessage("e-values not found."));

  return e_values_re;
}



Real DenseMatrix::determinant() const
{
  Assert(this->size1()==this->size2(), ExcMessage("Should be square"));

  using namespace boost::numeric::ublas;
  using PMatrix = permutation_matrix<int>;

  const int n = this->size1();
  BoostMatrix A(*this);
  PMatrix P(n);

  int res = lu_factorize(A, P);
  AssertThrow(res == 0, ExcMessage("LU factorization failed!"));

  Real det = 1.;
  int sign = 1;
  for (int i=0; i<n; ++i)
  {
    det *= A(i,i);
    if (P(i) != i)
      sign *= -1;
  }

  return sign * det;
}


//todo should be const
DenseMatrix
DenseMatrix::
inverse(Real &det) const
{
  using namespace boost::numeric::ublas;
  using namespace boost::numeric::ublas;
  using PMatrix = permutation_matrix<int>;

  const int n = this->size1();
  BoostMatrix A(*this);
  PMatrix P(n);

  int res = lu_factorize(A,P);

  AssertThrow(res == 0, ExcMessage("LU factorization failed!"));

  // create identity matrix of "inverse"
  BoostMatrix inv_A = identity_matrix<Real>(n);

  // backsubstitute to get the inverse
  lu_substitute(A, P, inv_A);

  det = 1.;
  int sign = 1;
  for (int i=0; i<n; ++i)
  {
    det *= A(i,i);
    if (P(i) != i)
      sign *= -1;
  }

  return inv_A;
}



Size
DenseMatrix::
get_num_rows() const
{
  return this->size1();
}

Size
DenseMatrix::
get_num_cols() const
{
  return this->size2();
}

Real
DenseMatrix::
norm_frobenius() const
{
  const Size n_rows = this->get_num_rows();
  const Size n_cols = this->get_num_cols();

  Real norm = 0.0;
  for (Index row = 0 ; row < n_rows ; ++row)
    for (Index col = 0 ; col < n_cols ; ++col)
    {
      const Real value = (*this)(row,col);
      norm += value * value;
    }

  return sqrt(norm);
}

Real
DenseMatrix::
norm_max() const
{
  const Size n_rows = this->get_num_rows();
  const Size n_cols = this->get_num_cols();

  Real norm = 0.0;
  for (Index row = 0 ; row < n_rows ; ++row)
    for (Index col = 0 ; col < n_cols ; ++col)
    {
      const Real value = (*this)(row,col);
      norm = std::max(std::abs(value),norm);
    }

  return norm;
}

Real
DenseMatrix::
norm_infinity() const
{
  const Size n_rows = this->get_num_rows();
  const Size n_cols = this->get_num_cols();

  Real norm = 0.0;
  for (Index row = 0 ; row < n_rows ; ++row)
  {
    Real sum = 0.0;
    for (Index col = 0 ; col < n_cols ; ++col)
    {
      const Real value = (*this)(row,col);
      sum += std::abs(value);
    }
    norm = std::max(sum,norm);
  }
  return norm;
}

Real
DenseMatrix::
norm_one() const
{
  const Size n_rows = this->get_num_rows();
  const Size n_cols = this->get_num_cols();

  Real norm = 0.0;
  for (Index col = 0 ; col < n_cols ; ++col)
  {
    Real sum = 0.0;
    for (Index row = 0 ; row < n_rows ; ++row)
    {
      const Real value = (*this)(row,col);
      sum += std::abs(value);
    }
    norm = std::max(sum,norm);
  }
  return norm;
}

bool
DenseMatrix::
is_symmetric() const
{
  bool is_symmetric = true;

  const Size n_rows = this->get_num_rows();
  const Size n_cols = this->get_num_cols();
  for (int row = 0 ; row < n_rows ; ++row)
  {
    for (int col = row+1 ; col < n_cols ; ++col)
    {
      if ((*this)(row,col) != (*this)(col,row))
      {
        is_symmetric = false;
        break;
      }
    }

  }

  return is_symmetric;
}


void
DenseMatrix::
print_info(LogStream &out) const
{
  const auto size1 = this->size1();
  const auto size2 = this->size2();

  out << '[' << size1 << ',' << size2 << "](";
  if (size1 > 0)
  {
    out << '(' ;
    if (size2 > 0)
      out << (*this)(0, 0);
    for (int j = 1; j < size2; ++ j)
      out << ',' << (*this)(0, j);
    out << ')';
  }
  for (int i = 1; i < size1; ++ i)
  {
    out << ",(" ;
    if (size2 > 0)
      out << (*this)(i, 0);
    for (int j = 1; j < size2; ++ j)
      out << ',' << (*this)(i, j);
    out << ')';
  }
  out << ')';
}



void eig_dense_matrix(const DenseMatrix &A,
                      SafeSTLVector<Real> &eigenvalues_real,
                      SafeSTLVector<Real> &eigenvalues_imag,
                      DenseMatrix &eigenvectors)
{
  const int n_rows = A.get_num_rows();
#ifndef NDEBUG
  const int n_cols = A.get_num_cols();
  Assert(n_rows == n_cols,ExcDimensionMismatch(n_rows,n_cols));


  Assert(eigenvalues_real.size() == n_rows,ExcDimensionMismatch(eigenvalues_real.size(),n_rows));
  Assert(eigenvalues_imag.size() == n_rows,ExcDimensionMismatch(eigenvalues_imag.size(),n_rows));

  Assert(eigenvectors.get_num_rows() == n_rows,ExcDimensionMismatch(eigenvectors.get_num_rows(),n_rows));
  Assert(eigenvectors.get_num_cols() == n_cols,ExcDimensionMismatch(eigenvectors.get_num_cols(),n_cols));
#endif
  const int n = n_rows;
  DenseMatrix eigenvectors_trans(n,n);

  Teuchos::LAPACK<int, double> lapack;
  int info;

  const int workspace_size = 4*n; // this should be >= 3*n (higher values gives better performances)
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



  if (std::is_sorted(eigenvalues_real.begin(),eigenvalues_real.end()))
  {
    eigenvectors = boost::numeric::ublas::trans(eigenvectors_trans);
  }
  else
  {
    std::map<Real,std::set<Index>> eigval_real_sort;
    for (int i = 0 ; i < n ; ++i)
      eigval_real_sort[eigenvalues_real[i]].insert(i);

    const auto eigenvalues_imag_tmp = eigenvalues_imag;
    int i = 0;
    for (const auto &uniqueval_ids : eigval_real_sort)
    {
      const Real value = uniqueval_ids.first;
      for (const auto id : uniqueval_ids.second)
      {
        eigenvalues_real[i] = value;
        eigenvalues_imag[i] = eigenvalues_imag_tmp[id];

        // transposing and reordering the eigenvectors
        for (int row = 0 ; row < n ; ++row)
          eigenvectors(row,i) = eigenvectors_trans(id,row);

        ++i;
      } // end loop ids with same eigenvalue
    }
  }
}


void eig_dense_matrix_symm(const DenseMatrix &A,
                           SafeSTLVector<Real> &eigenvalues,
                           DenseMatrix &eigenvectors)
{
  const int n_rows = A.get_num_rows();
#ifndef NDEBUG
  const int n_cols = A.get_num_cols();
  Assert(n_rows == n_cols,ExcDimensionMismatch(n_rows,n_cols));
  Assert(A.is_symmetric(),ExcMessage("The matrix is not symmetric."));


  Assert(eigenvalues.size() == n_rows,ExcDimensionMismatch(eigenvalues.size(),n_rows));

  Assert(eigenvectors.get_num_rows() == n_rows,ExcDimensionMismatch(eigenvectors.get_num_rows(),n_rows));
  Assert(eigenvectors.get_num_cols() == n_cols,ExcDimensionMismatch(eigenvectors.get_num_cols(),n_cols));
#endif

  Teuchos::LAPACK<int, double> lapack;
  int info;

  const int n = n_rows;
  const int workspace_size = 10*n; // thsi should be >= 3*n (higher values gives better performances)
  double workspace[workspace_size];

  // The Lapack SYEV routine assumes that the matrix entries are sorted column-wise
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


  if (std::is_sorted(eigenvalues.begin(),eigenvalues.end()))
  {
    eigenvectors = boost::numeric::ublas::trans(eigenvectors_trans);
  }
  else
  {
    // sorting the eigenvalues in ascending order
    std::map<Real,std::set<Index>> eigval_real_sort;
    for (int i = 0 ; i < n ; ++i)
      eigval_real_sort[eigenvalues[i]].insert(i);

    int i = 0;
    for (const auto &uniqueval_ids : eigval_real_sort)
    {
      const Real value = uniqueval_ids.first;
      for (const auto id : uniqueval_ids.second)
      {
        eigenvalues[i] = value;

        // transposing and reordering the eigenvectors
        for (int row = 0 ; row < n ; ++row)
          eigenvectors(row,i) = eigenvectors_trans(id,row);

        ++i;
      } // end loop ids with same eigenvalue
    }
  }
}

#if 0
#ifdef SERIALIZATION

template<class Archive>
void
DenseMatrix::
serialize(Archive &ar, const unsigned int version)
{
  ar &boost::serialization::make_nvp("DenseMatrix_base_t",
                                     boost::serialization::base_object<BoostMatrix>(*this));
}
#endif //SERIALIZATION
#endif
IGA_NAMESPACE_CLOSE


#if 0
#ifdef SERIALIZATION

BOOST_CLASS_EXPORT_IMPLEMENT(iga::DenseMatrix)
template void iga::DenseMatrix::serialize(OArchive &, const unsigned int);
template void iga::DenseMatrix::serialize(IArchive &, const unsigned int);

#endif // SERIALIZATION
#endif
