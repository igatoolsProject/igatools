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

#ifndef __IGA_ARRAY_UTILS_H_
#define __IGA_ARRAY_UTILS_H_

#include <igatools/base/config.h>
#include <igatools/utils/safe_stl_array.h>

#include <array>
#include <algorithm>
#include <utility>

IGA_NAMESPACE_OPEN



/**
 * Returns an array filled with the sequence of <tt>N</tt> integers numbers
 * from <tt>init</tt> to <tt>init+N-1</tt>.
 */
template<int N>
constexpr
auto
sequence(const int init = 0)
-> SafeSTLArray<int, N>
{
    SafeSTLArray<int, N> seq;
    std::iota(seq.begin(),seq.end(),init);
    return seq;
}


IGA_NAMESPACE_CLOSE

#endif
