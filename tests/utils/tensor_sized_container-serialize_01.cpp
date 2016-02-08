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

/**
 *  @file
 *  @brief  TensorSizedContainer serialization
 *  @author  martinelli
 *  @date 2015-05-05
 */

#include "../tests.h"
#include <igatools/utils/tensor_sized_container.h>


template <int N>
void tsc_serialization()
{
  OUTSTART

  using T = TensorSizedContainer<N>;
//  using T = Test<2>;

  T size(3);

  const std::string filename = "tsc.xml";

  {
    ofstream os(filename);
    OArchive archive(os);
    archive << size;
  }
  out.begin_item("TensorSizedContainer<" + std::to_string(N) + "> before serialization");
  size.print_info(out);
  out.end_item();

  T size_tmp;
  {
    ifstream is(filename);
    IArchive archive(is);
    archive >> size_tmp;
  }
  out.begin_item("TensorSizedContainer<" + std::to_string(N) + "> after serialization");
  size_tmp.print_info(out);
  out.end_item();

  OUTEND
}



int main()
{
//  tsc_serialization<0>();
  tsc_serialization<1>();
  tsc_serialization<2>();
  tsc_serialization<3>();
  return 0;
}
