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
 *  Test for BSplineVecSpace, constructing a component space
 *
 *  author: pauletti
 *  date: 2015-03-17
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>



template<int dim, int codim=0>
void component_space(const int deg=3,  const int n_knots = 10)
{
  OUTSTART

  auto grid  = CartesianGrid<dim>::create(n_knots);
  using VecSpace = BSplineSpace<dim, dim+codim>;
  using CompSpace  = BSplineSpace<dim, 1>;
  typename VecSpace::Degrees degt(deg);
  typename VecSpace::Periodicity periodic(false);
  periodic[0] = true;
  typename VecSpace::EndBehaviour end_b(BasisEndBehaviour::interpolatory);
  end_b[0] = BasisEndBehaviour::periodic;

  auto space = VecSpace::create(degt, grid, InteriorReg::maximum, periodic, end_b);
  space->print_info(out);

  auto comp_space = CompSpace::create(space->get_degree_table()[0],grid, InteriorReg::maximum,
                                      space->get_periodicity()[0],
                                      space->get_end_behaviour_table()[0]);

  comp_space->print_info(out);

  OUTEND
}




int main()
{

  component_space<2, 1>();

  return 0;
}
