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

/**
 *  @file
 *  @brief  SafeSTLVector serialization
 *  @author  martinelli
 *  @date 2015-05-05
 */

#include "../tests.h"
#include <igatools/utils/safe_stl_vector.h>

void vector_serialization()
{
  OUTSTART

  SafeSTLVector<Real> vec = {1.,2.,3.,4.,5.};


  {
    ofstream os("vector.xml");
    OArchive archive(os);
    archive << vec;
  }
  out.begin_item("SafeSTLVector<Real> before serialization");
  vec.print_info(out);
  out.end_item();


  vec.clear();
  {
    ifstream is("vector.xml");
    IArchive archive(is);
    archive >> vec;
  }
  out.begin_item("SafeSTLVector<Real> after serialization");
  vec.print_info(out);
  out.end_item();

  OUTEND
}



int main()
{
  vector_serialization();
  return 0;
}
