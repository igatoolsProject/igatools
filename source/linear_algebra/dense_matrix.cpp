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

IGA_NAMESPACE_CLOSE



#ifdef SERIALIZATION
/*
 * The next macro is needed to instantiate the serialize() function for the active archives
 */
BOOST_CLASS_EXPORT(iga::DenseMatrix)

#endif // SERIALIZATION
