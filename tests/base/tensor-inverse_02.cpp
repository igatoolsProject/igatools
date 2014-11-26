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
    auto B = inverse(A, det);
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
    compute_inverse<4,4>();

    compute_inverse<1,2>();
    compute_inverse<2,1>();

    compute_inverse<3,2>();
    compute_inverse<2,3>();

    return 0;
}
