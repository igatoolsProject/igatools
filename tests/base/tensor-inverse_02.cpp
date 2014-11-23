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
//TODO: Add standard description

/*
 *  Developing a new function for handling tensor inverses
 *
 *  author: pauletti
 *  date: 2014-11-23
 *
 */

#include "../tests.h"
#include <igatools/base/tensor.h>

template<typename T>
using Inverse =
    Conditional<
    T::is_tensor,
    Conditional<
    T::value_t::rank==0,
    T,
    Tensor<T::value_t::dim,T::value_t::rank,typename T::value_t::tensor_t::co_type,
    Tensor<T::dim,T::rank,typename T::tensor_t::co_type,typename T::value_t::value_t> >
    >,
    Tdouble >;


template<class T>
inline
EnableIf<(T::dim==0) && (T::value_t::dim == 0) &&
         (T::rank==1) && (T::value_t::rank == 1), Inverse<T> >
new_inverse_(const T &A, Real &det)
{
	det = 1.;
    return Inverse<T>();
}

template<class T>
inline
EnableIf<(T::dim==1) && (T::value_t::dim == 1) &&
         (T::rank==1) && (T::value_t::rank == 1), Inverse<T> >
new_inverse_(const T &A, Real &det)
{
	det = A[0][0];
   	Assert(det != Real(0.0), ExcDivideByZero());

   	Inverse<T> A_inv;
    A_inv[0][0] =  1.0 / det;

    return A_inv;
}



template<class T>
inline
EnableIf<(T::dim==2) && (T::value_t::dim == 2) &&
         (T::rank==1) && (T::value_t::rank == 1), Inverse<T> >
new_inverse_(const T &A, Real &det)
{
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    Assert(det != Real(0.0), ExcDivideByZero());

    Inverse<T> A_inv;
    const Real InvDet = 1.0 / det;

    A_inv[0][0] = A[1][1] * InvDet;
    A_inv[0][1] = A[0][1] * (-InvDet);
    A_inv[1][0] = A[1][0] * (-InvDet);
    A_inv[1][1] = A[0][0] * InvDet;

    return A_inv;
}



template<class T>
inline
EnableIf<(T::dim==3) && (T::value_t::dim == 3) &&
         (T::rank==1) && (T::value_t::rank == 1), Inverse<T> >
new_inverse_(const T &A, Real &det)
{
	Inverse<T> A_inv;

    const Real t4 = A[0][0]*A[1][1];
    const Real t6 = A[0][0]*A[1][2];
    const Real t8 = A[0][1]*A[1][0];
    const Real t00 = A[0][2]*A[1][0];
    const Real t01 = A[0][1]*A[2][0];
    const Real t04 = A[0][2]*A[2][0];
    det = (t4*A[2][2]-t6*A[2][1]-t8*A[2][2]+
    		t00*A[2][1]+t01*A[1][2]-t04*A[1][1]);
    Assert(det != Real(0.0), ExcDivideByZero());

    const Real t07 = 1.0/det;
    A_inv[0][0] = (A[1][1]*A[2][2]-A[1][2]*A[2][1])*t07;
    A_inv[0][1] = (A[0][2]*A[2][1]-A[0][1]*A[2][2])*t07;
    A_inv[0][2] = (A[0][1]*A[1][2]-A[0][2]*A[1][1])*t07;
    A_inv[1][0] = (A[1][2]*A[2][0]-A[1][0]*A[2][2])*t07;
    A_inv[1][1] = (A[0][0]*A[2][2]-t04)*t07;
    A_inv[1][2] = (t00-t6)*t07;
    A_inv[2][0] = (A[1][0]*A[2][1]-A[1][1]*A[2][0])*t07;
    A_inv[2][1] = (t01-A[0][0]*A[2][1])*t07;
    A_inv[2][2] = (t4-t8)*t07;

    return  A_inv;
}



template<class T>
EnableIf<(T::dim == T::value_t::dim), Inverse<T> >
new_inverse(const T &A, Real &det)
{
	return new_inverse_(A, det);
}


/**
 * right inverse, has the property that:
 * A*R = I for all a in the range of A and
 * R*A = I
 */
template<class T>
EnableIf<(T::dim < T::value_t::dim), Inverse<T> >
new_inverse(const T &A, Real &det)
{
	const auto A_t   = co_tensor(transpose(A));
    const auto G     = compose(A_t, A);
    const auto G_inv = new_inverse(G, det);
	det = sqrt(det);

	return compose(G_inv, A_t);
}



template<class T>
EnableIf<(T::dim > T::value_t::dim), Inverse<T> >
new_inverse(const T &A, Real &det)
{
	const auto A_t   = co_tensor(transpose(A));
    const auto G     = compose(A, A_t);
    const auto G_inv = new_inverse(G, det);
	det = sqrt(det);

	return compose(A_t, G_inv);
}


template<int rdim, int cdim>
void compute_inverse()
{
	OUTSTART

    Tensor<cdim, 1, tensor::covariant, Tensor< rdim, 1, tensor::contravariant, Tdouble> > A;

    for (int i = 0; i < cdim; ++i)
    	for (int j = 0; j < rdim; ++j)
    		A[i][j] = cos(i*j);

    out << "A =" << endl;
    out << A << endl;

    Real det;
    auto B = new_inverse(A, det);
    out << "Determinant(A) = " << det << endl;

    out << "A^(-1) =" << endl;
    out << B << endl;

    out << "A o A^(-1) =" << endl;
    out << compose(A,B)  << endl;
    out << "A^(-1) o A =" << endl;
    out << compose(B,A)  << endl;

    OUTEND
}


int main()
{
	compute_inverse<1,1>();
	compute_inverse<2,2>();
	compute_inverse<3,3>();

	compute_inverse<1,2>();
	compute_inverse<2,1>();

	compute_inverse<3,2>();
	compute_inverse<2,3>();

    return 0;
}
