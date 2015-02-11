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


#ifndef TOPOLOGY_INLINE_H_
#define TOPOLOGY_INLINE_H_

#include <igatools/geometry/topology.h>
#include <igatools/geometry/unit_element.h>

IGA_NAMESPACE_OPEN


template<int dim>
inline
TopologyId<dim>::
TopologyId(const Index id)
    :
    id_(id)
{
    Assert(id_>=-1,ExcLowerRange(id_,-1));
}

template<int dim>
inline
Index
TopologyId<dim>::
get_id() const
{
    return id_;
}

template<int dim>
inline
bool
TopologyId<dim>::
is_element() const
{
    return (id_ == -1)?true:false;
}

template<int dim>
inline
bool
TopologyId<dim>::
is_face() const
{
    return (id_ >= 0)?true:false;
}


template<int dim>
inline
vector<Index>
TopologyId<dim>::
get_active_directions() const
{
    vector<Index> active_directions;
    if (this->is_element())
    {
        for (Index d = 0 ; d < dim ; ++d)
            active_directions.emplace_back(d);

        Assert(active_directions.size() == dim,ExcDimensionMismatch(active_directions.size(),dim));
    }
    else if (this->is_face())
    {
        const int face_id = this->get_id();
        for (const Index &d : UnitElement<dim>::face_active_directions[face_id])
            active_directions.emplace_back(d);

        Assert(active_directions.size() == UnitElement<dim>::face_dim,
               ExcDimensionMismatch(active_directions.size(),UnitElement<dim>::face_dim));
    }
    else
    {
        Assert(false,ExcInvalidState());
        AssertThrow(false,ExcInvalidState());
    }

    return active_directions;
}



template<int dim>
inline
ElemTopology<dim>::
ElemTopology()
    :
    TopologyId<dim>(-1)
{}


template<int dim>
inline
FaceTopology<dim>::
FaceTopology(const Index face_id)
    :
    TopologyId<dim>(face_id)
{
    Assert(face_id >= 0 && face_id < UnitElement<dim>::faces_per_element,
           ExcIndexRange(face_id,0,UnitElement<dim>::faces_per_element));
};



IGA_NAMESPACE_CLOSE




#endif // TOPOLOGY_INLINE_H_
