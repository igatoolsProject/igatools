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
/*
 * Test for the determinant of square tensors.
 * cavallini 2012/11/7
 */

#include "../tests.h"

#include <igatools/base/tensor.h>

template< int dim >
void det_eval()
{
    Tensor<dim, 1, tensor::covariant, Tensor <dim,1, tensor::contravariant, Tdouble > > t;

    double accumulator = 0;

    out << "The determinant of:" << std::endl;

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim ; j++)
        {
            accumulator = accumulator + 1;
            t[i][j] = accumulator;

            out << "t["<< i << "]" << "["<< j << "] = "
                << t[i][j] << std::endl;
        }
    }

    double det;

    det = determinant<dim,dim>(t);

    out << "is " << det << std::endl;

}

int main(int argc, char *argv[])
{

    det_eval<1>();
    det_eval<2>();
    det_eval<3>();

}

