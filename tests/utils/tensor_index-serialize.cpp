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
 *  @brief  TensorIndex serialization
 *  @author  martinelli
 *  @date 2015-09-08
 */

#include "../tests.h"
#include <igatools/utils/tensor_index.h>

/*
template <int M>
class Test
{
public:
  Test() = default;
  Test(const int N)
  :
    N_(N)
  {
    array_.fill(N);
  }

private:
  int N_;

  SafeSTLArray<int,M> array_;

  friend class cereal::access;

  template<class Archive>
  void serialize(Archive &ar)
  {
    ar & make_nvp("N_",N_);
    ar & make_nvp("array_",array_);
  }
};
//*/

template <int N>
void ti_serialization()
{
  OUTSTART

  using T = TensorIndex<N>;

  T index(3);

  const std::string filename = "ti.xml";

  {
    ofstream os(filename);
    OArchive archive(os);
    archive << index;
  }
  out.begin_item("TensorIndex<" + std::to_string(N) + "> before serialization");
  index.print_info(out);
  out.end_item();

  T index_tmp;
  {
    ifstream is(filename);
    IArchive archive(is);
    archive >> index_tmp;
  }
  out.begin_item("TensorIndex<" + std::to_string(N) + "> after serialization");
  index_tmp.print_info(out);
  out.end_item();

  OUTEND
}



int main()
{
  ti_serialization<0>();
  ti_serialization<1>();
  ti_serialization<2>();
  ti_serialization<3>();
  return 0;
}
