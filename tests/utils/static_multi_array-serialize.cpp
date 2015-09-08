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
 * Test for StaticMultiArray: serialization
 * author: martinelli
 * data:   08 Sep 2015
 *
 */

#include "../tests.h"
#include <igatools/utils/static_multi_array.h>

template <int dim,int rank = 2>
void static_multiarray_serialization()
{
  OUTSTART

  using T = StaticMultiArray<int,dim,rank>;

  T data(3);
  data.fill_progression();

  using std::to_string;
  const std::string name = "StaticMultiArray<int," + to_string(dim) + "," + to_string(rank) + ">";

  out.begin_item(name + " before serialization");
  data.print_info(out);
  out.end_item();


  const std::string filename = "static_multiarray_" +
                               to_string(dim) + "_" + to_string(rank) + ".xml";

  {
    ofstream os(filename);
    OArchive archive(os);
    archive << data;
  }

  T data_tmp;
  {
    ifstream is(filename);
    IArchive archive(is);
    archive >> data_tmp;
  }
  out.begin_item(name + " after serialization");
  data_tmp.print_info(out);
  out.end_item();


  OUTEND
}


int main()
{
  out.depth_console(10);

  static_multiarray_serialization<1>();
  static_multiarray_serialization<2>();
//  static_multiarray_serialization<3>();

  return 0;
}
