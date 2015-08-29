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


// QualityAssurance: martinelli, 04 Sep 2014

#include <igatools/utils/unique_id_generator.h>
#include <igatools/base/exceptions.h>

#include <limits>


IGA_NAMESPACE_OPEN


Index
UniqueIdGenerator::
get_unique_id()
{
  const Index int_max = std::numeric_limits<Index>::max();
  Assert(id_ < int_max,ExcIndexRange(id_,0,int_max-1));
  AssertThrow(id_ < int_max,ExcIndexRange(id_,0,int_max-1));

  Index ret = id_;
  id_++;
  return ret;
}


Index UniqueIdGenerator::id_ = 0;

IGA_NAMESPACE_CLOSE

