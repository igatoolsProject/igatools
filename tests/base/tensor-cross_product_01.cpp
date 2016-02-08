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
 *  cross product
 *
 *  author: pauletti
 *  date: 2014-11-23
 *
 */

#include "../tests.h"
#include <igatools/base/tensor.h>


template <int dim, int codim = 1>
void
compute_cp()
{
  OUTSTART
  Derivatives<dim, dim+1, 1, 1> DF;
  for (int i = 0; i < dim; ++i)
  {
    DF[i][i] = 1;
  }

  auto cp = cross_product<dim, codim>(DF);
  out << "cross product: " << cp << endl;

  for (int i = 0; i < dim; ++i)
  {
    out << "scalar product: " << i << ": " << scalar_product(cp, DF[i]) << endl;
  }
  OUTEND
}



template <int dim, int codim = 1>
void
compute_cp2()
{
  OUTSTART
  Derivatives<dim, dim+1, 1, 1> DF;
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim+1; ++j)
      DF[i][j] = i+j;
  }

  auto cp = cross_product<dim, codim>(DF);
  out << "cross product: " << cp << endl;

  for (int i = 0; i < dim; ++i)
  {
    out << "scalar product: " << i << ": " << scalar_product(cp, DF[i]) << endl;
  }
  OUTEND
}


int main()
{
  compute_cp<1>();
  compute_cp<2>();

  compute_cp2<1>();
  compute_cp2<2>();

  return 0;
}
