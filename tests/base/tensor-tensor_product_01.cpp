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
 * Test for the
 * Author: pauletti 2013/04/26
 *
 */

#include "../tests.h"

#include <igatools/base/tensor.h>


template<int dim, int rdim>
void test_tensor_product()
{

    Points<dim>  a;
    Points<rdim> b;

    out << "The tensor product of:" << std::endl;
    out << a << endl;
    out << "and:" << endl;
    out << b << endl;
    out << "is:" << endl;
    out << tensor_product(a, b) << endl;
}

int main(int argc, char *argv[])
{
    out.depth_console(1);

    test_tensor_product<2,2>();
    test_tensor_product<2,3>();

    return 0;
}

