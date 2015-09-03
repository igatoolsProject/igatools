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
 *  Test for the SplineSpace class  subspace extraction facilities.
 *  get_sub_space_mult and get_sub_space_degree
 *
 *  author: pauletti
 *  date: 2014-11-07
 */

#include "../tests.h"
#include <igatools/geometry/unit_element.h>
#include <igatools/basis_functions/spline_space.h>


template<int k, int dim, int range=1, int rank=1>
void sub_space(const TensorSize<dim> &n_knots, const TensorIndex<dim> &degree)
{
  OUTSTART
  using SplineSpace = SplineSpace<dim, range, rank>;
  auto grid = Grid<dim>::create(n_knots);
  typename SplineSpace::DegreeTable deg {degree};
  auto int_mult = SplineSpace::get_multiplicity_from_regularity(InteriorReg::maximum,
                  deg, grid->get_num_intervals());
  auto space = SplineSpace::create(deg, grid, int_mult);

  for (auto  s_id : UnitElement<dim>::template elems_ids<k>())
  {
    out.begin_item("Sub element id: " + std::to_string(s_id));

    out.begin_item("Multiplicity");
    auto sub_mult = space->template get_sub_space_mult<k>(s_id);
    sub_mult.print_info(out);
    out.end_item();

    out.begin_item("Degree");
    auto sub_deg = space->template get_sub_space_degree<k>(s_id);
    sub_deg.print_info(out);
    out.end_item();

    out.begin_item("Periodicity");
    auto sub_periodic = space->template get_sub_space_periodicity<k>(s_id);
    //sub_periodic.print_info(out);
    out.end_item();

    out.end_item();
  }
  OUTEND
}



int main()
{
  out.depth_console(10);

  sub_space<1,2>(TensorSize<2>({3,4}), TensorIndex<2>({1,3}));

  return  0;
}
