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

#include <igatools/basis_functions/function_space.h>
#include <igatools/utils/unique_id_generator.h>

IGA_NAMESPACE_OPEN

template <class Grid_>
constexpr  std::array<Size, Grid_::dim> FunctionSpaceOnGrid<Grid_>::dims;

template <class Grid_>
FunctionSpaceOnGrid<Grid_>::
FunctionSpaceOnGrid(std::shared_ptr<GridType> grid)
    :
    GridWrapper<GridType>(grid),
    id_(UniqueIdGenerator::get_unique_id())
{};


template <class GridType>
Index
FunctionSpaceOnGrid<GridType>::
get_id() const
{
    return id_;
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/function_space.inst>
