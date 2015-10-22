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

#ifndef TYPES_H_
#define TYPES_H_

#include <igatools/base/config.h>

IGA_NAMESPACE_OPEN

/**
* Bit field flags for specifying which VTK grid type must be used.
*/
enum class vtkGridType : std::int64_t
{
    /** VTK structured grid */
    Structured          =    1 << 0,

    /** VTK unstructure grid with linear cells */
    UnstructuredLinear  =    1 << 1,

    /** VTK unstructure grid with quadratic cells */
    UnstructuredQuadratic  = 1 << 2,

    /** None */
    None                   = 1 << 3,
};

IGA_NAMESPACE_CLOSE

#endif // TYPES_H_
