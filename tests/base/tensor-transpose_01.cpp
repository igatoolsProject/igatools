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
/*
 *  Test for computing the transpose of a tensor.
 *
 *  author: pauletti
 *  date: Feb 18, 2013
 *
 */

#include "../tests.h"

#include <igatools/base/tensor.h>


template <int dim1, int rank1, class type1, int dim2, int rank2, class type2>
void run_test1()
{

    out << "Test for the transpose of a tensor:" << endl;
    out << "on a tensor of ranks: " << rank1 << ", "<< rank2 << endl;
    out << " and dims: " << dim1 <<", " << dim2 << ".\n";

    typedef Tensor< dim1, rank1, type1, Tensor<dim2,rank2,type2,Tdouble> > T_t;

    T_t A;

    const int size1 = T_t::size;
    const int size2 = T_t::value_t::size;

    for (int i=0; i<size1; ++i)
        for (int j=0; j<size2; ++j)
        {
            A(i)(j) = size2*i + j;
        }



    auto B = transpose(A);
    out << "A= "<< endl << A << endl;
    out << "A^t= "<< endl << B << endl;


    out << endl;
}


int main()
{
    out.depth_console(10);

    run_test1<1, 1, tensor::covariant, 2, 1, tensor::contravariant>();
    run_test1<2, 1, tensor::covariant, 2, 1, tensor::contravariant>();
    run_test1<2, 1, tensor::covariant, 2, 2, tensor::contravariant>();

    run_test1<3, 1, tensor::covariant, 3, 1, tensor::contravariant>();
    run_test1<3, 1, tensor::covariant, 2, 1, tensor::contravariant>();



    return  0;
}

