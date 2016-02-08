//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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
 * Test for contract 1
 * Author: pauletti 2013/04/26
 *
 */
#include "../tests.h"

#include <igatools/base/tensor.h>




template<int dim, int rdim>
void test_contract_1()
{
  const int order = 2;
  Derivatives <dim, rdim, 1, order> D2F;
  Derivatives <dim, rdim, 1, order-1> DF;

  double val = 0;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < rdim; ++j)
      DF[i][i] = ++val;

  Real det;
  auto DF_inv = inverse(DF, det);


  out << "The contract_1 of:" << endl;
  out << D2F << endl;
  out << "with:" << endl;
  out << DF << endl;
  out << "is:" << endl;
  out << contract_1(D2F,DF) << endl;

  out << "The contract_1 of:" << endl;
  out << D2F << endl;
  out << "with:" << endl;
  out << co_tensor(transpose(DF_inv)) << endl;
  out << "is:" << endl;
  out << contract_1(D2F,co_tensor(transpose(DF_inv))) << endl;
  out << endl;
}




int main()
{

  out.depth_console(1);

  test_contract_1<1,1>();
  test_contract_1<1,3>();
  test_contract_1<1,3>();

  test_contract_1<2,2>();
  test_contract_1<2,3>();

  test_contract_1<3,3>();


  return 0;

}

