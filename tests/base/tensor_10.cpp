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
 *  Test for redesign tensor class.
 *  Here test for the basic features like constructor and index access.
 *
 *  author: pauletti
 *  date: Feb 15, 2013
 *
 */

#include "../tests.h"

#include <igatools/base/tensor.h>



//Define some tensors
template <int dim, int rank, class type>
void run_test1()
{
    out << "Test for standard Constructor" << endl;
    out << "Tensor of rank: " << rank << " and dim: " << dim << endl;

    out<< "A =" << endl;
    Tensor< dim, rank, type, Tdouble > A;
    out << A << endl;

    out<< "B =" << endl;
    Tensor< dim, rank, typename type::co_type, Tdouble > B;
    out << B << endl;

    out<< "C =" << endl;
    Tensor< dim, rank, type, Tensor< dim, rank, typename type::co_type, Tdouble > > C;
    out << C << endl;

    out<< "D =" << endl;
    Tensor< dim, rank, typename type::co_type, Tensor< dim, rank, type, Tdouble > > D;
    out << D << endl;
    out << endl;

    out << "A+A=" << endl;
    out << A+A << endl;
}


//Test linear index access
template <int dim, class type>
void run_test2()
{
    const int rank = 2;
    out << "Test for filling a tensor" << endl;
    out << "Tensor of rank: " << rank << " " << dim << endl;

    Tensor< dim, rank, type, Tdouble > A;

    for (int i=0; i<dim; ++i)
        for (int j=0; j<dim; ++j)
            A[i][j] = i*dim + j;

    out << A << endl;
    out << endl;
}

//Test linear index access
template <int dim, class type>
void run_test3()
{
    const int rank = 2;
    out << "Test accesing unsing tensor indices" << endl;
    out << "Tensor of rank: " << rank << " and dim: " << dim << endl;

    Tensor< dim, rank, type, Tdouble > A;


    for (int i=0; i<dim; ++i)
        for (int j=0; j<dim; ++j)
            A[i][j] = i*dim + j;
    out << A << endl;


    TensorIndex<rank> index;
    for (int i=0; i<dim; ++i)
        for (int j=0; j<dim; ++j)
        {
            index[0] = i;
            index[1] = j;
            out << A(index) << "  ";
            out << A[i][j]  << endl;
        }
    out << endl;
}


int main()
{
    out.depth_console(10);

//  run_test1<3, 1, tensor::covariant>();

    run_test1<1, 1, tensor::covariant>();
    run_test1<3, 1, tensor::covariant>();

    run_test1<1,2, tensor::covariant>();
    run_test2<1, tensor::covariant>();
    run_test3<1, tensor::covariant>();

    run_test1<2,2, tensor::covariant>();
    run_test2<2, tensor::covariant>();
    run_test3<2, tensor::covariant>();

    return  0;
}
