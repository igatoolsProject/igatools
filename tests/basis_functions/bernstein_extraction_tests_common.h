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

/**
 *  @file
 *  @brief Common code for tests/basis_functions/bernstein_extraction_*.cpp
 *  @author martinelli
 *  @date Feb 10, 2016
 */

#include "../tests.h"
#include <igatools/basis_functions/spline_space.h>
#include <igatools/basis_functions/space_tools.h>


template <int dim,int range>
typename SplineSpace<dim,range>::
template ComponentContainer<SafeSTLArray<SafeSTLVector<Real>,dim>>
compute_knots_with_repetitions(
  const SplineSpace<dim,range> &spline_space,
  typename SplineSpace<dim,range>::BoundaryKnotsTable bdry_knots,
  typename SplineSpace<dim,range>::EndBehaviourTable ebt)
{
  using SplineSpace = SplineSpace<dim,range>;


  const auto &grid = *spline_space.get_grid();
  const auto &degree_table = spline_space.get_degree_table();
  const auto &multiplicity_table = spline_space.get_interior_mult();
  const auto &periodicity_table = spline_space.get_periodicity();

  const auto active_comp_ids = spline_space.get_active_components_id();

  using Knots = SafeSTLArray<SafeSTLVector<Real>,dim>;
  using KnotsTable = typename SplineSpace::template ComponentContainer<Knots>;
  KnotsTable knots_with_repetitions;

  for (int comp : active_comp_ids)
  {
    const auto &degree_table_comp = degree_table[comp];

    auto &knots_with_repetitions_comp = knots_with_repetitions[comp];
    knots_with_repetitions_comp =
      space_tools::compute_knots_with_repetition(
        grid,
        degree_table_comp,
        multiplicity_table[comp],
        ebt[comp],
        bdry_knots[comp],
        periodicity_table[comp]);
  }

  return knots_with_repetitions;
}


