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
 *  Test for the BallFunction class as a mapping
 *  Computing the volume
 *
 *  author: pauletti
 *  date: 2014-10-24
 */

#include "../tests.h"

#include <igatools/functions/function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/functions/function_lib.h>
#include <igatools/functions/identity_function.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>


template <int dim>
Real ball_volume(const int n_knots)
{
  using Function = functions::BallFunction<dim>;

  auto flag = ValueFlags::w_measure|ValueFlags::point;

  auto quad = QGauss<dim>(3);


  BBox<dim> box;
  box[0] = {0., 1.};
  for (int i=1; i<dim-1; ++i)
    box[i] = {0., M_PI};
  if (dim>1)
    box[dim-1] = {0., 2. * M_PI};

  auto grid = Grid<dim>::create(box, n_knots);
  auto F = Function::create(grid, IdentityFunction<dim>::create(grid));

  using Mapping   = Mapping<dim, 0>;
  using ElementIt = typename Mapping::ElementIterator;
  auto map = Mapping::create(F);
  map->template reset<dim>(flag, quad);

//    ElementIt elem(map, 0);
//    ElementIt end(map, IteratorState::pass_the_end);
  ElementIt elem = map->begin();
  ElementIt end = map->end();

  map->template init_cache<dim>(elem);
  Real vol = 0.;
  for (; elem != end; ++elem)
  {
    map->template fill_cache<dim>(elem, 0);
    const auto w_meas = elem->template get_w_measures<dim>(0);

    for (auto &w : w_meas)
      vol += w;
  }
  return vol;
}

template<int dim>
void
comp_volume()
{
  OUTSTART
  for (int n_knots = 2; n_knots < 6; ++n_knots)
    out << n_knots << "\t" << ball_volume<dim>(n_knots) << endl;

  auto exact = std::pow(M_PI, dim/2.) / tgamma(dim/ 2. + 1);
  out << "Exact volume: " <<  exact << endl;
  OUTEND
}

int main()
{
  out.depth_console(10);

  comp_volume<1>();
  comp_volume<2>();
  comp_volume<3>();

  return 0;
}
