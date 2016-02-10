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

#ifndef SPACE_TOOLS_H_
#define SPACE_TOOLS_H_

#include <igatools/geometry/grid.h>

IGA_NAMESPACE_OPEN

/**
 * @brief This namespace contains some useful functions for operations involving the
 * data that are used to define a SplineSpace.
 */
namespace space_tools
{


SafeSTLVector<Real>
compute_knots_with_repetition_1D(
  const SafeSTLVector<Real> &knots_no_repetitions,
  const int deg_comp_dir,
  const SafeSTLVector<int> &mult_comp_dir,
  const BasisEndBehaviour &ends_comp_dir,
  const SafeSTLArray<SafeSTLVector<Real>,2> &boundary_knots_comp_dir,
  const bool &periodic_comp_dir)
{
#ifndef NDEBUG
  const auto &l_knots = boundary_knots_comp_dir[0];
  const auto &r_knots = boundary_knots_comp_dir[1];

  if (periodic_comp_dir)
  {
    Assert(ends_comp_dir == BasisEndBehaviour::periodic,
           ExcMessage("Periodic inconsistency"));
    Assert(l_knots.size()==0,
           ExcMessage("Periodic inconsistency"));
    Assert(r_knots.size()==0,
           ExcMessage("Periodic inconsistency"));
  }
  else
  {
    if (ends_comp_dir == BasisEndBehaviour::interpolatory)
    {
      Assert(l_knots.size()==0,
             ExcMessage("Interpolatory inconsistency"));
      Assert(r_knots.size()==0,
             ExcMessage("Interpolatory inconsistency"));
    }
    if (ends_comp_dir == BasisEndBehaviour::end_knots)
    {
//      const auto &knots = grid_->get_knot_coordinates(j);
      const int m = deg_comp_dir + 1;
      Assert(l_knots.size() == m,
             ExcMessage("Wrong number of boundary knots"));
      Assert(r_knots.size() == m,
             ExcMessage("Wrong number of boundary knots"));
      Assert(knots_no_repetitions.front() >= l_knots.back(),
             ExcMessage("Boundary knots should be smaller or equal a"));
      Assert(knots_no_repetitions.back() <= r_knots.front(),
             ExcMessage("Boundary knots should be greater or equal b"));
      Assert(std::is_sorted(l_knots.begin(), l_knots.end()),
             ExcMessage("Boundary knots is not sorted"));
      Assert(std::is_sorted(r_knots.begin(), r_knots.end()),
             ExcMessage("Boundary knots is not sorted"));
    }
  }
#endif


  const auto order = deg_comp_dir + 1;

//    const int m = order;
  const int K = std::accumulate(mult_comp_dir.begin(),mult_comp_dir.end(),0);


  SafeSTLVector<Real> knots_comp_dir(2*order+K);

  auto rep_it = knots_comp_dir.begin() + order;
  auto m_it = mult_comp_dir.begin();
  auto k_it = ++knots_no_repetitions.begin();
  auto end = mult_comp_dir.end();
  for (; m_it !=end; ++m_it, ++k_it)
  {
    for (int iMult = 0; iMult < *m_it; ++iMult, ++rep_it)
      *rep_it = *k_it;
  }


  if (periodic_comp_dir)
  {
    const Real a = knots_no_repetitions.front();
    const Real b = knots_no_repetitions.back();
    const Real L = b - a;
    for (int i = 0 ; i < order ; ++i)
    {
      knots_comp_dir[i] = knots_comp_dir[K+i] - L;
      knots_comp_dir[K+order+i] = knots_comp_dir[order+i] + L;
    }
  }
  else
  {
    const auto &endb = ends_comp_dir;
    switch (endb)
    {
      case BasisEndBehaviour::interpolatory:
      {
        const Real a = knots_no_repetitions.front();
        const Real b = knots_no_repetitions.back();

        for (int i = 0 ; i < order ; ++i)
        {
          knots_comp_dir[i] = a;
          knots_comp_dir[order+K+i] = b;
        }
      }
      break;
      case BasisEndBehaviour::end_knots:
      {
        const auto &left_knts = boundary_knots_comp_dir[0];
        const auto &right_knts = boundary_knots_comp_dir[1];

        for (int i = 0 ; i < order ; ++i)
        {
          knots_comp_dir[i] = left_knts[i];
          knots_comp_dir[K+order+i] = right_knts[i];
        }
      }
      break;
      case BasisEndBehaviour::periodic:
        Assert(false, ExcMessage("Impossible"));
    }
  }

  return knots_comp_dir;
}



template<int dim_>
SafeSTLArray<SafeSTLVector<Real>,dim_>
compute_knots_with_repetition(
  const Grid<dim_> &grid,
  const SafeSTLArray<int,dim_> &deg_comp,
  const SafeSTLArray<SafeSTLVector<int>,dim_> &mult_comp,
  const SafeSTLArray<BasisEndBehaviour,dim_> &ends_comp,
  const SafeSTLArray<SafeSTLArray<SafeSTLVector<Real>,2>,dim_> &boundary_knots_comp,
  const SafeSTLArray<bool,dim_> &periodic_comp)
{
  SafeSTLArray<SafeSTLVector<Real>,dim_> knots_comp;
  for (int dir = 0 ; dir < dim_ ; ++dir)
  {
    knots_comp[dir] = compute_knots_with_repetition_1D(
                        grid.get_knot_coordinates(dir),
                        deg_comp[dir],
                        mult_comp[dir],
                        ends_comp[dir],
                        boundary_knots_comp[dir],
                        periodic_comp[dir]);
  }

  return knots_comp;
}


}

IGA_NAMESPACE_CLOSE

#endif
