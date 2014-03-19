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


#ifndef TOPOLOGY_INLINE_H_
#define TOPOLOGY_INLINE_H_

#include <igatools/geometry/topology.h>


IGA_NAMESPACE_OPEN


inline
TopologyId::
TopologyId(const Index id)
    :
    id_(id)
{
    Assert(id_>=-1,ExcLowerRange(id_,-1));
}

inline
Index
TopologyId::
get_id() const
{
    return id_;
}

inline
bool
TopologyId::
is_element() const
{
    return (id_ == -1)?true:false;
}

inline
bool
TopologyId::
is_face() const
{
    return (id_ >= 0)?true:false;
}


inline
ElemTopology::
ElemTopology()
    :
    TopologyId(-1)
{}


inline
FaceTopology::
FaceTopology(const Index id)
    :
    TopologyId(id)
{
    Assert(id>=0,ExcMessage("Face ID must be positive."));
};



IGA_NAMESPACE_CLOSE




#endif // TOPOLOGY_INLINE_H_
