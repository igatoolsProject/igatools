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
 *  Test for the l2_projection function.
 *  Bspline spaces case
 *
 *  author: pauletti
 *  date: 2013-10-10
 */

#include "../tests.h"
#include "./common_functions.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/formula_function.h>

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>

#include <igatools/basis_functions/space_tools.h>


template<int dim , int range ,int rank, LAPack la_pack>
void project_l2(const int p, const int num_knots = 10)
{
  OUTSTART

  using Space =  BSplineSpace<dim,range,rank>;
  using RefSpace =  ReferenceSpace<dim,range,rank>;

  auto knots = Grid<dim>::create(num_knots);
  auto space = Space::create(p, knots) ;

  const int n_qpoints = 4;
  QGauss<dim> quad(n_qpoints);

  auto f = BoundaryFunction<dim>::create(knots);
  auto proj_func = space_tools::projection_l2<RefSpace,la_pack>(f, space, quad);
  proj_func->print_info(out);

  OUTEND
}



int main()
{
#if defined(USE_TRILINOS)
  const auto la_pack = LAPack::trilinos_epetra;
#elif defined(USE_PETSC)
  const auto la_pack = LAPack::petsc;
#endif

  project_l2<0,1,1, la_pack>(1);
  project_l2<1,1,1, la_pack>(3);
  project_l2<2,1,1, la_pack>(3);
  project_l2<3,1,1, la_pack>(1);

  return 0;
}

