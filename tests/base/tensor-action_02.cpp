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
 * Test for the action  of rank 2 tensors
 * Author: pauletti 2013/04/26
 *
 */

#include "../tests.h"

#include <igatools/base/tensor.h>


template<int dim, int rdim>
void test_action()
{
    const int order = 2;
    Derivatives <dim, rdim, 1, order> D2F;
    Points<dim> u;
//
//  for (int i = 0; i < cdim; i++)
//  {
//      t2[i] = i;
//      for (int j = 0; j < rdim; j++)
//          t1[i][j] = cdim*i+j;
//  }

    out << "The action of:" << std::endl;
    out << D2F << endl;
    out << "on:" << endl;
    out << u << std::endl;
    out << "is:" << std::endl;
    out << action(D2F, u) << endl;
}

int main(int argc, char *argv[])
{
    out.depth_console(1);
//    do_action<1,1>();
//  do_action<2,3>();
    test_action<2,2>();
//  do_action<3,2>();
//  do_action<2,3>();
//  do_action<3,3>();

}

