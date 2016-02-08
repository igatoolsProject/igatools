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
 *  Test for refinement of the Grid
 *
 *  author: martinelli
 *  date: 2013
 *
 */

#include "../tests.h"

#include <igatools/geometry/grid.h>



template <int dim>
void test_evaluate()
{
  using Grid = Grid<dim>;
  auto grid = Grid::create();
//    const std::string active_property = "active";
  const std::string influence_property = "influence";

//    grid->add_property(active_property);
  grid->add_property(influence_property);


  typename Grid::IndexType elem_id(0,TensorIndex<dim>());
//    grid->get_element_property(active_property).insert(elem_id);
  grid->get_elements_with_property(influence_property).emplace_back(elem_id);


  out << "===============================================================" << endl;
  out << "D i m = " << dim << endl;

  out << "---------------------------------------------------------------" << endl;
  out << "Unrefined Grid:" << endl;
  grid->print_info(out);
  out << "---------------------------------------------------------------" << endl;
  out << endl;

  out << "---------------------------------------------------------------" << endl;
  out << "Previous grid refined with 2 sub_intervals along each direction:" << endl;
  grid->refine();
  grid->print_info(out);
  out << "---------------------------------------------------------------" << endl;
  out << endl;

  out << "---------------------------------------------------------------" << endl;
  out << "Previous grid refined with 3 sub_intervals along each direction:" << endl;
  grid->refine(3);
  grid->print_info(out);
  out << "---------------------------------------------------------------" << endl;
  out << endl;

  out << "---------------------------------------------------------------" << endl;
  out << "Previous grid refined with 4 sub_intervals along the direction 0:" << endl;
  grid->refine_direction(0,4);
  grid->print_info(out);
  out << "---------------------------------------------------------------" << endl;
  out << "===============================================================" << endl;

  out << endl;

  out << endl;

}

int main()
{
  out.depth_console(10);

  test_evaluate<1>();
  out << endl ;

  test_evaluate<2>();
  out << endl ;

  test_evaluate<3>();
  out << endl ;

  return 0;
}
