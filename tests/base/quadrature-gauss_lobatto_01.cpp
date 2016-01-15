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
 *  Test for  Gauss-Lobatto quadrature scheme
 *
 *  author: pauletti
 *  date: 2015-03-06
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/array_utils.h>

template <int dim>
void uniform(const int n_pts)
{
  OUTSTART

  QGaussLobatto<dim> quad(n_pts);
  quad.print_info(out) ;

  OUTEND
}
template <int dim>
void non_uniform(const int n_pts)
{
  OUTSTART

  const TensorSize<dim> n_points(sequence<dim>(n_pts));
  QGaussLobatto<dim> quad(n_points);
  quad.print_info(out);

  OUTEND
}


int main()
{
  const int max_n_pts = 5;
  for (int n = 2; n <  max_n_pts; ++n)
  {
    uniform<0>(n);
    uniform<1>(n);
    uniform<2>(n);
    uniform<3>(n);

    non_uniform<0>(n);
    non_uniform<1>(n);
    non_uniform<2>(n);
    non_uniform<3>(n);
  }

  return 0;
}



