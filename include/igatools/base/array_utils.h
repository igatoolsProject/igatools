//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
#include <array>
#include <algorithm>

IGA_NAMESPACE_OPEN

template <class T, int dim>
inline
std::array<T, dim> filled_array(const T &v)
{
    std::array<T,dim> res;
    res.fill(v);
    return res;
}



/**
 * Returns an array filled with the sequence of numbers
 */
template <int dim, class T = int>
inline
std::array<T,dim>
sequence(const int init = 0)
{
    std::array<T,dim> res;
    std::iota(res.begin(), res.end(), init);
    return res;
}

IGA_NAMESPACE_CLOSE

#endif
