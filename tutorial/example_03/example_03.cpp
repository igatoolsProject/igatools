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

#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_element.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
// [quad include]
#include <igatools/base/quadrature_lib.h>
// [quad include]
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

template <int dim>
void loop_on_grid_with_cache()
{
  // [loop as before]
  out << "Traversing the elements of a " << dim << "-dimensional grid." << endl;
  const int n_knots = 3;
  auto grid = Grid<dim>::const_create(n_knots);

  auto cache_handler = grid->create_cache_handler();

  auto quad = QGauss<dim>::create(2);
  auto flag = grid_element::Flags::weight;

  cache_handler->set_element_flags(flag);

  auto elem = grid->begin();
  const auto elem_end = grid->end();
  // [loop as before]
  // [init cache]
  cache_handler->init_element_cache(elem,quad);
  // [init cache]

  int elem_id = 0;
  for (; elem != elem_end; ++elem)
  {
    // [fill cache]
    cache_handler->fill_element_cache(elem);
    // [fill cache]
    out << "The tensor index of element: " << elem_id << " is: "<< elem->get_index() << endl;

    // [get meas]
    auto w_meas = elem->template get_weights<dim>(0);
    out << "The weighted measure is: ";
    w_meas.print_info(out);
    // [get meas]
    out << endl;

    ++elem_id;
  }
  out << endl;
}


template <int dim>
void loop_on_space_with_cache()
{
  out << "Traversing the elements of a " << dim << "-dimensional B-spline space." << endl;
  const int n_knots = 3;
  auto grid = Grid<dim>::const_create(n_knots);
  const int degree = 2;
  auto space = SplineSpace<dim>::const_create(degree, grid);
  auto basis = BSpline<dim>::const_create(space);

  auto cache_handler = basis->create_cache_handler();
  auto quad = QGauss<dim>::create(1);
  auto flag = space_element::Flags::value;

  cache_handler->set_element_flags(flag);

  auto elem = basis->begin();
  const auto elem_end = basis->end();
  cache_handler->init_element_cache(elem,quad);

  for (; elem != elem_end; ++elem)
  {
    cache_handler->fill_element_cache(elem);
    out << "Element: " << elem->get_index() << " has global basis: ";
    elem->get_local_to_global(DofProperties::active).print_info(out);
    out << endl;

    out.begin_item("Basis values:");
    elem->template get_basis_element<space_element::_Value>().print_info(out);
    out.end_item();
  }
  out << endl;
}


int main()
{

  loop_on_grid_with_cache<1>();
  loop_on_grid_with_cache<2>();
  loop_on_grid_with_cache<3>();

  loop_on_space_with_cache<1>();
  loop_on_space_with_cache<2>();
  loop_on_space_with_cache<3>();

  return 0;
}



