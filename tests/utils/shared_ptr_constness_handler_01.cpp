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
 * Test for SafeSTLVector, SafeSTLArray and their serialization
 * author: pauletti, martinelli
 * date:   2014-08-26
 * date:   2015-05-05 (added the tests for serialization)
 *
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/utils/shared_ptr_constness_handler.h>





using namespace iga;

using std::shared_ptr;

using Grid = Grid<1>;

const std::string class_name = "SharedPtrConstnessHandler";

void do_test_nonconst()
{
  OUTSTART

  shared_ptr<Grid> grid_nonconst = Grid::create();

  SharedPtrConstnessHandler<Grid> constness_handler(grid_nonconst);

  out.begin_item(class_name + "::print_info()");
  constness_handler.print_info(out);
  out.end_item();

  out.begin_item(class_name + "::get_ptr_data()");
  constness_handler.get_ptr_data()->print_info(out);
  out.end_item();

  out.begin_item(class_name + "::get_ptr_const_data()");
  constness_handler.get_ptr_const_data()->print_info(out);
  out.end_item();

  out.begin_item(class_name + "::get_ref_data()");
  constness_handler.get_ref_data().print_info(out);
  out.end_item();

  out.begin_item(class_name + "::get_ref_const_data()");
  constness_handler.get_ref_const_data().print_info(out);
  out.end_item();

#if 0
  out.begin_item(class_name + "::operator*()");
  Grid &grid_ref = *constness_handler;
  grid_ref.print_info(out);
  out.end_item();

  out.begin_item(class_name + "::operator->()");
  constness_handler->print_info(out);
  out.end_item();
#endif

  OUTEND
}


void do_test_const()
{
  OUTSTART

  shared_ptr<const Grid> grid_const = Grid::create();

  SharedPtrConstnessHandler<Grid> constness_handler(grid_const);

  out.begin_item(class_name + "::print_info()");
  constness_handler.print_info(out);
  out.end_item();

  out.begin_item(class_name + "::get_ptr_const_data()");
  constness_handler.get_ptr_const_data()->print_info(out);
  out.end_item();

  out.begin_item(class_name + "::get_ref_const_data()");
  constness_handler.get_ref_const_data().print_info(out);
  out.end_item();

  out.begin_item(class_name + "::operator*()");
  (*constness_handler).print_info(out);
  out.end_item();

  out.begin_item(class_name + "::operator->()");
  constness_handler->print_info(out);
  out.end_item();

  OUTEND
}


int main()
{
  OUTSTART

  do_test_nonconst();

  do_test_const();

  OUTEND

  return 0;

}
