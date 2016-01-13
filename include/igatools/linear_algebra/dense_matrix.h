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


#ifndef DENSE_MATRIX_H_
#define DENSE_MATRIX_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/safe_stl_vector.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

IGA_NAMESPACE_OPEN

/**
 * @brief Dense matrix with real entries.
 *
 * A typical use of it is the local matrix.
 * @note This class inherits (also the constructors) from the matrix type in the
 * <a href="http://www.boost.org/doc/libs/1_55_0/libs/numeric/ublas/doc/index.htm">Boost uBlas library</a>
 * Therefore it can be used as a boost::numeric::ublas::matrix<Real>.
 *
 * @ingroup linear_algebra
 * @ingroup serializable
 *
 * @author M. Martinelli, 2012, 2013, 2014
 * @author S. Pauletti, 2012, 2013, 2014
 */
class DenseMatrix : public boost::numeric::ublas::matrix<Real>
{
public:
  /** Type of the base class. */
  using BoostMatrix =  boost::numeric::ublas::matrix<Real> ;

  /** We inherith the constructors of the base class. */
  using BoostMatrix::BoostMatrix;
#if 0
  DenseMatrix(const Index n, const bool init_to_zero = true);
  DenseMatrix(const Index n_rows, const Index n_cols, const bool init_to_zero = true);

  DenseMatrix(const DenseMatrix &matrix) = default;
  DenseMatrix(DenseMatrix &&matrix) = default;
#endif

  /**
   * Assignment operator for assigning zeros to all entries of the matrix
   * by writing
   * @code
     DenseMatrix matrix;
     ... // working with the matrix
     matrix = 0.0; // reset the matrix entries to zero
     @endcode
   * @note If used in Debug mode with a @p value different from zero,
   * an assertion will be raised.
   */
  DenseMatrix &operator=(const Real value) ;

  /** Type for a row of the matrix. */
  using MatrixRowType = boost::numeric::ublas::matrix_row<const BoostMatrix> ;

  /** Returns the mtrix row identified by @p row_id. */
  MatrixRowType
  get_row(const int row_id) const;

  /** Returns the number of rows in the DenseMatrix. */
  int size1() const;

  /** Returns the number of columns in the DenseMatrix. */
  int size2() const;


  /**
   * Matrix inversion routine.
   * @note Uses lu_factorize and lu_substitute in uBLAS to invert a matrix
   */
  DenseMatrix inverse(Real &det) const;

  SafeSTLVector<Real> eigen_values() const;

  Real determinant() const;

  /**
   * Returns the <em>Frobenius norm</em> of the matrix, i.e.
   * \f[
     || A ||_F = \biggl( \sum_{i,j} |a_{ij}|^2 \biggr)^{1/2}
     \f]
   * where \f$ a_{ij} \f$ is the
   * entry in the \f$i\f$-th row and \f$j\f$-th column of the matrix \f$ A\f$.
   */
  Real norm_frobenius() const;

  /**
   * Returns the <em>max-norm</em> (or <em>uniform norm</em>) of the matrix, i.e.
   * \f[
     || A ||_{\text{max}} = \max_{i,j} |a_{ij}|
     \f]
   * where \f$ a_{ij} \f$ is the
   * entry in the \f$i\f$-th row and \f$j\f$-th column of the matrix \f$ A\f$.
   */
  Real norm_max() const;

  /**
   * Returns the <em>one-norm</em> of the matrix, i.e.
   * \f[
     || A ||_1 = \max_{j} \sum_{i} |a_{ij}|
     \f]
   * where \f$ a_{ij} \f$ is the
   * entry in the \f$i\f$-th row and \f$j\f$-th column of the matrix \f$ A\f$.
   */
  Real norm_one() const;

  /**
   * Returns the <em>infinity-norm</em> of the matrix, i.e.
   * \f[
     || A ||_{\infty} = \max_{i} \sum_{j} |a_{ij}|
     \f]
   * where \f$ a_{ij} \f$ is the
   * entry in the \f$i\f$-th row and \f$j\f$-th column of the matrix \f$ A\f$.
   */
  Real norm_infinity() const;

  /** Returns the number of rows of the matrix. */
  Size get_num_rows() const;

  /** Returns the number of columns of the matrix. */
  Size get_num_cols() const;

  void print_info(LogStream &out) const;

  /**
   * Returns true if the matrix entries verify \f$ (A)_{ij} = (A)_{ji} \quad \forall i \, , \, j = 1,\dots,n\f$
   * @return
   */
  bool is_symmetric() const;

  bool operator==(const DenseMatrix &in) const
  {
    bool is_equal = (this == &in);
    if (!is_equal)
    {
      bool equal_size = (this->get_num_rows() == in.get_num_rows() &&
                         this->get_num_cols() == in.get_num_cols());
      if (equal_size)
      {
        const auto &this_data_begin = this->data().begin();
        const auto &this_data_end   = this->data().end();
        const auto &in_data_begin = in.data().begin();
        is_equal = std::equal(this_data_begin,this_data_end,in_data_begin);
      }
      else
      {
        is_equal = false;
      }
    }
    return is_equal;
  }

private:

#ifdef SERIALIZATION
  /**
   * @name Functions needed for serialization
   */
  ///@{
  friend class cereal::access;

  template <class Archive>
  void
  save(Archive &ar) const;

  template <class Archive>
  void
  load(Archive &ar);
  ///@}
#endif // SERIALIZATION
};



/**
 * Computes the eigenvalues (and the associated eigenvectors) of the matrix @p A.
 *
 * The eigenvalues are sorted in ascending order relatively to their real part and for the i-th eigenvalue,
 * the associated eigenvector is the i-th column of the matrix @p eigenvectors.
 *
 * @relates DenseMatrix
 */
void eig_dense_matrix(const DenseMatrix &A,
                      SafeSTLVector<Real> &eigenvalues_real,
                      SafeSTLVector<Real> &eigenvalues_imag,
                      DenseMatrix &eigenvectors);

/**
 * Computes the eigenvalues (and the associated eigenvectors) of the symmetric matrix @p A.
 *
 * The @p eigenvalues are sorted in ascending order and for the i-th eigenvalue, the associated eigenvector is
 * the i-th column of the matrix @p eigenvectors.
 *
 * @note The input matrix A must be symmetric, otherwise an assertion will be raised (in Debug mode).
 *
 * @relates DenseMatrix
 */
void eig_dense_matrix_symm(const DenseMatrix &A,
                           SafeSTLVector<Real> &eigenvalues,
                           DenseMatrix &eigenvectors);


IGA_NAMESPACE_CLOSE


#ifdef SERIALIZATION

CEREAL_SPECIALIZE_FOR_ARCHIVE(IArchive,iga::DenseMatrix,cereal::specialization::member_load_save)
CEREAL_SPECIALIZE_FOR_ARCHIVE(OArchive,iga::DenseMatrix,cereal::specialization::member_load_save)

#endif // SERIALIZATION

#endif /* DENSE_MATRIX_H_ */
