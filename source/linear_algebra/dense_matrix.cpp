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

IGA_NAMESPACE_OPEN

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

IGA_NAMESPACE_CLOSE

