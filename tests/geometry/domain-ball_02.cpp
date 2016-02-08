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
 *  Test for the BallGridFunction class as a domain
 *  Computing the volume
 *
 *  author: pauletti
 *  date: 2014-10-24
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/geometry/domain_handler.h>


template <int dim>
Real ball_volume(const int n_knots)
{
  BBox<dim> box;
  box[0] = {0., 1.};
  for (int i=1; i<dim-1; ++i)
    box[i] = {0., M_PI};
  if (dim>1)
    box[dim-1] = {0., 2. * M_PI};

  auto grid = Grid<dim>::const_create(box, n_knots);

  using Ball = grid_functions::BallGridFunction<dim>;
  auto ball = Ball::const_create(grid);
  auto domain = Domain<dim,0>::const_create(ball);

  auto domain_handler = domain->create_cache_handler();

  auto elem = domain->begin();
  auto end = domain->end();


  using Flags = domain_element::Flags;
  auto flag = Flags::w_measure;

  domain_handler->set_element_flags(flag);

  auto quad = QGauss<dim>::create(3);

  domain_handler->init_element_cache(elem,quad);
  Real vol = 0.;
  for (; elem != end; ++elem)
  {
    domain_handler->fill_element_cache(elem);
    const auto w_meas = elem->get_element_w_measures();

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
