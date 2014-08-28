//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
#include <boost/numeric/ublas/lu.hpp>

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


DenseMatrix
DenseMatrix::
inverse()
{
    using namespace boost::numeric::ublas;
    using pmatrix_t = permutation_matrix<std::size_t>;

    // create a working copy of the input
    BoostMatrix A(static_cast<BoostMatrix &>(*this)) ;

    // create a permutation matrix for the LU-factorization
    pmatrix_t P(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A,P);

    AssertThrow(res == 0, ExcMessage("LU factorization failed!"));

    // create identity matrix of "inverse"
    BoostMatrix inv_A = identity_matrix<Real>(A.size1());

    // backsubstitute to get the inverse
    lu_substitute(A, P, inv_A);

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


IGA_NAMESPACE_CLOSE

