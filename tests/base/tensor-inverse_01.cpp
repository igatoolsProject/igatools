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
// Test for the inverse of a rank 2 Tensor and other things too
// Author: Seba 2012/10/26

#include "../tests.h"

#include <igatools/base/tensor.h>


template<int rdim, int cdim>
void do_inverse()
{

    Tensor<cdim, 1, tensor::covariant, Tensor< rdim, 1, tensor::contravariant, Tdouble> > A;
    Tensor<rdim, 1, tensor::covariant, Tensor< cdim, 1, tensor::contravariant, Tdouble> > B;
    Tensor<cdim, 1, tensor::covariant, Tensor< rdim, 1, tensor::contravariant, Tdouble> > C;

    Tensor< cdim, 1, tensor::contravariant, Tdouble> u;
    Tensor< rdim, 1, tensor::contravariant, Tdouble> v;

    for (int i = 0; i < cdim; i++)
    {
        u[i] = i+1;
        for (int j = 0; j < rdim; j++)
            A[i][j] = cos(i*j)   ;
    }

    out << "Case:" << rdim << " " << cdim << std::endl;
    out << "The Tensor A is:" << std::endl;
    out << A << std::endl;

    out << "Its action on:" << std::endl;
    out << u << std::endl;
    out << "is:" << std::endl;
    v = action(A,u);
    out << v << std::endl;


//    out << "Its transpose action on:" << std::endl;
//    out << v << std::endl;
//    out << "is:" << std::endl;
//    u = transpose_action(A,v);
//    out << u << std::endl;




    out << "Its inverse is:" << std::endl;
    inverse<cdim, rdim> (A,B);
    out << B << std::endl;

    out << "compose A by its inverse to get:" << std::endl;
    out << compose(A,B)  << std::endl;

    out << std::endl;
}

int main(int argc, char *argv[])
{
    do_inverse<2,2>();
    do_inverse<3,2>();
    do_inverse<3,3>();
    //do_inverse<2,1>();

}

