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

#include <igatools/basis_functions/bspline.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/nurbs.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/functions/ig_function.h>
#include <igatools/io/writer.h>

// [includes]
#include <igatools/geometry/grid_element.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_handler.h>
#include <igatools/base/quadrature_lib.h>
// [includes]

using namespace iga;
using namespace std;

LogStream out;

// [quarter_annulus]
shared_ptr<const Domain<2>> quarter_annulus(const Size nel)
{
  using numbers::PI;
  BBox<2> box;
  box[0] = {{1.0,2.0}};
  box[1] = {{0.0,PI/2}};
  auto grid       = Grid<2>::const_create(box,nel+1);
  auto geom_funct = grid_functions::BallGridFunction<2>::const_create(grid);
  return Domain<2>::const_create(geom_funct);
}
// [quarter_annulus]

// [grid_loop]
template <int dim>
void grid_loop(shared_ptr<const Grid<dim>> grid)
{

  auto grid_el     = grid->begin();
  auto grid_el_end = grid->end();

  out << "Traversing the elements of a " << dim << "-dimensional grid." << endl;
  for (; grid_el!=grid_el_end; ++grid_el)
  {
    out << "The flat/tensor indices of the current element are:  " << grid_el->get_index() << endl;
  }
  out << endl;
}
// [grid_loop]

// [basis_loop_start]
template <int dim, int codim=0, int range=1, int rank=1>
void basis_loop_with_cache(shared_ptr<const Basis<dim,codim,range,rank>> basis)
{

  auto basis_el      = basis->begin();
  auto basis_el_end  = basis->end();
  auto cache_handler = basis->create_cache_handler();
// [basis_loop_start]

// [basis_loop_set]
  auto flag = basis_element::Flags::value;
  cache_handler->set_element_flags(flag);
// [basis_loop_set]

// [basis_loop_init]
  auto quad = QGauss<dim>::create(1);
  cache_handler->init_element_cache(basis_el,quad);
// [basis_loop_init]

// [basis_loop_loop]
  for (; basis_el!=basis_el_end; ++basis_el)
  {
    cache_handler->fill_element_cache(basis_el);
// [basis_loop_loop]

// [basis_loop_view]
    auto vals = basis_el->get_element_values();

    out.begin_item("Basis values:");
    vals.print_info(out);
    out.end_item();
  }
// [basis_loop_view]
}

// [main_trivial]
int main()
{
  auto segment = Grid<1>::const_create(3);
  grid_loop<1>(segment);

  auto square = Grid<2>::const_create(3);
  grid_loop<2>(square);

  auto cube = Grid<3>::const_create(3);
  grid_loop<3>(cube);
// [main_trivial]

// [main_basis_loop_ref]
  auto space = SplineSpace<2>::const_create(2,square);
  auto basis = BSpline<2>::const_create(space);

  out << "Traversing basis functions on the reference domain:" << endl;
  basis_loop_with_cache<2>(dynamic_pointer_cast<const Basis<2,0,1,1>>(basis));
// [main_basis_loop_ref]

// [main_basis_loop_phy]
  auto annulus   = quarter_annulus(2);
  auto grid      = annulus->get_grid_function()->get_grid();
  auto spl_space = SplineSpace<2>::const_create(2,grid);
  auto ref_basis = BSpline<2>::const_create(spl_space);
  auto phy_basis = PhysicalBasis<2>::const_create(ref_basis,annulus);

  out << "Traversing basis functions on the annulus domain:" << endl;
  basis_loop_with_cache<2>(dynamic_pointer_cast<const Basis<2,0,1,1>>(phy_basis));
// [main_basis_loop_phy]

  return 0;
}
