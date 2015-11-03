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

// [func_grid]
template <int dim>
void loop_on_grid_with_cache()
{
  // [func_grid]
  // [loop_as_before]
  out << "Traversing the elements of a " << dim << "-dimensional grid." << endl;
  const int n_knots = 3;
  auto grid = Grid<dim>::const_create(n_knots);

  auto elem = grid->begin();
  const auto elem_end = grid->end();
  // [loop_as_before]

  // [create_handler]
  auto cache_handler = grid->create_cache_handler();
  // [create_handler]

  // [set_cache]
  const auto flag = grid_element::Flags::weight;
  cache_handler->set_element_flags(flag);
  // [set_cache]

  // [init_cache]
  auto quad = QGauss<dim>::create(2);
  cache_handler->init_element_cache(elem,quad);
  // [init_cache]

  int elem_id = 0;
  for (; elem != elem_end; ++elem)
  {
    // [fill_cache]
    cache_handler->fill_element_cache(elem);
    // [fill_cache]

    out << "The tensor index of element: " << elem_id << " is: "<< elem->get_index() << endl;

    // [get_meas]
    const auto &w_meas = elem->get_element_weights();
    out << "The weighted measure is: ";
    w_meas.print_info(out);
    // [get_meas]
    out << endl;

    ++elem_id;
  }
  out << endl;
}


// [func_basis]
template <int dim>
void loop_on_basis_with_cache()
{
  // [func_basis]
  out << "Traversing the elements of a " << dim << "-dimensional B-spline space." << endl;
  const int n_knots = 3;
  auto grid = Grid<dim>::const_create(n_knots);
  const int degree = 2;
  auto space = SplineSpace<dim>::const_create(degree, grid);
  auto basis = BSpline<dim>::const_create(space);

  auto elem_basis = basis->begin();
  const auto elem_basis_end = basis->end();


  auto cache_handler = basis->create_cache_handler();

  // [set_basis_flags]
  auto flag = space_element::Flags::value;
  cache_handler->set_element_flags(flag);
  // [set_basis_flags]


  auto quad = QGauss<dim>::create(1);
  cache_handler->init_element_cache(elem_basis,quad);

  for (; elem_basis != elem_basis_end; ++elem_basis)
  {
    cache_handler->fill_element_cache(elem_basis);

    out << "Element: " << elem_basis->get_index() << " has global basis: ";
    elem_basis->get_local_to_global().print_info(out);
    out << endl;

    // [basis_values]
    out.begin_item("Basis values:");
    elem_basis->get_element_values().print_info(out);
    out.end_item();
    // [basis_values]
  }
  out << endl;
}


int main()
{

  loop_on_grid_with_cache<1>();
  loop_on_grid_with_cache<2>();
  loop_on_grid_with_cache<3>();

  loop_on_basis_with_cache<1>();
  loop_on_basis_with_cache<2>();
  loop_on_basis_with_cache<3>();

  return 0;
}



