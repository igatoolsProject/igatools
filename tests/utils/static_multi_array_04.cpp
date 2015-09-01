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
 * Test for StaticMultiArray
 * author: pauletti
 * data:   28 Aug 2014
 *
 */

#include "../tests.h"
#include <igatools/utils/static_multi_array.h>

template <int dim>
void fill_progression()
{
  OUTSTART

  StaticMultiArray<Index,dim,2> data1(4);
  out << "Fill progression from 2 "<< endl;
  data1.fill_progression(2);
  out << "Size: " << data1.flat_size() << endl;
  data1.print_info(out);
  out << endl;

  OUTEND
}


int main()
{
  out.depth_console(10);


  fill_progression<1>();
  fill_progression<2>();
  fill_progression<3>();

  return 0;
}
