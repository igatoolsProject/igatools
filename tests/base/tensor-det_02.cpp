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
 * Test for the determinant of square tensors.
 * cavallini 2012/11/7
 * updated by pauletti 2013/05/07
 */

#include "../tests.h"

#include <igatools/base/tensor.h>

template< int dim_dom, int dim_range >
void det_eval()
{
  Derivatives<dim_dom, dim_range, 1, 1> t;

  out << "The determinant of:" << std::endl;

  double accumulator = 0;

  for (int i = 0; i < dim_dom; i++)
  {
    for (int j = 0; j < dim_range; j++)
    {
      t[i][j] = ++accumulator;
    }
  }

  out << t << endl;

  double det = determinant<dim_dom, dim_range>(t);

  out << "is " << det << std::endl;

}

int main(int argc, char *argv[])
{

  det_eval<1,2>();
  det_eval<1,3>();
  det_eval<2,3>();

}
