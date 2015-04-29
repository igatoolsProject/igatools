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


#include <igatools/geometry/grid_wrapper.h>

IGA_NAMESPACE_OPEN


template <class GridType>
GridWrapper<GridType>::
GridWrapper(std::shared_ptr<GridType> grid)
    :
    grid_ {grid}
{
    Assert(grid_ != nullptr,ExcNullPtr());
}


template <class GridType>
GridWrapper<GridType>::
~GridWrapper()
{
//    refine_h_connection_.disconnect();
    insert_knots_connection_.disconnect();
}


template <class GridType>
std::shared_ptr<GridType>
GridWrapper<GridType>::
get_grid()
{
    return grid_;
}


template <class GridType>
std::shared_ptr<const GridType>
GridWrapper<GridType>::
get_grid() const
{
    return grid_;
}


template <class GridType>
void
GridWrapper<GridType>::
refine_h_directions(
    const std::array<bool,Topology::dim> &refinement_directions,
    const std::array<Size,Topology::dim> &n_subdiv_directions)
{
    grid_->refine_directions(refinement_directions,n_subdiv_directions);
}


template <class GridType>
void
GridWrapper<GridType>::
refine_h_direction(const int direction_id, const Size n_subdivisions)
{
    grid_->refine_direction(direction_id,n_subdivisions);
}


template <class GridType>
void
GridWrapper<GridType>::
refine_h(const Size n_subdivisions)
{
    grid_->refine(n_subdivisions);
}

#if 0
template <class GridType>
void
GridWrapper<GridType>::
connect_refinement_h_function(const typename GridType::SignalRefineSlot &subscriber)
{
    refine_h_connection_ = grid_->connect_refinement(subscriber);
}
#endif

template <class GridType>
void
GridWrapper<GridType>::
connect_insert_knots_function(const typename GridType::SignalInsertKnotsSlot &subscriber)
{
    insert_knots_connection_ = grid_->connect_insert_knots(subscriber);
}


template <class GridType>
void
GridWrapper<GridType>::
insert_knots(SafeSTLArray<vector<Real>,GridType::dim> &knots_to_insert)
{
    grid_->insert_knots(knots_to_insert);
}


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_wrapper.inst>
