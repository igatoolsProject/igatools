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

#include <igatools/utils/tensor_range.h>

IGA_NAMESPACE_OPEN

template<>
SafeSTLSet<TensorIndex<1> >
el_tensor_range(TensorIndex<1> first, TensorIndex<1> last)
{
  Assert(first <= last, ExcMessage("first bigger than last"));
  SafeSTLSet<TensorIndex<1>> result;
  for (int i=first[0]; i<last[0]; ++i)
  {
    TensorIndex<1> el {i};
    result.insert(el);
  }

  return result;
}



template<>
SafeSTLSet<TensorIndex<0> >
el_tensor_range(TensorIndex<0> first, TensorIndex<0> last)
{
  SafeSTLSet<TensorIndex<0>> result;
  return result;
}

IGA_NAMESPACE_CLOSE
