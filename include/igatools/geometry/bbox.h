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

#ifndef BBOX_H_
#define BBOX_H_

#include <igatools/base/config.h>

#include <igatools/utils/safe_stl_array.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Bounding Box, a dim-dimensional rectangular
 * box described by the intervals.
 *
 * eg. BBox<2> box {{0,0.5},{1,2}}.
 */
template<int dim>
class BBox : public SafeSTLArray<SafeSTLArray<Real, 2>, dim>
{
public:
    using SafeSTLArray<SafeSTLArray<Real, 2>, dim>::SafeSTLArray;

    /**
     * Default constructor. It builds a box representing the <tt>dim</tt>-dimensional
     * unit hypercube \f$[0,1]^{dim}\f$.
     */
    BBox()
    {
        for (auto &bounds_dir : (*this))
        {
            bounds_dir[0] = 0.0;
            bounds_dir[1] = 1.0;
        }
    }
};


IGA_NAMESPACE_CLOSE


#endif // BBOX_H
