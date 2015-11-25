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

#include <igatools/base/objects_container.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/NURBS.h>

#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/container/map/map_fwd.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <boost/fusion/include/at_key.hpp>


using std::shared_ptr;
using namespace boost::fusion;

IGA_NAMESPACE_OPEN

ObjectsContainer::
ObjectsContainer()
{
};



auto
ObjectsContainer::
create() -> shared_ptr<self_t>
{
    return shared_ptr<self_t> (new self_t());
};



auto
ObjectsContainer::
const_create() -> shared_ptr<const self_t>
{
    return shared_ptr<const self_t> (new self_t());
};



template <class T>
void
ObjectsContainer::
insert_object (const std::shared_ptr<T> object,
               const Index &id)
{
    Assert (object != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert ((!this->is_object<T>(id)),
            ExcMessage("Object id already defined for the same object type"));

    auto &objects_T = at_key<T>(objects_);

    objects_T[id] = object;
};



template <class T>
auto
ObjectsContainer::
get_object (const Index &id) const -> std::shared_ptr<T>
{
    Assert ((this->is_object<T>(id)),
            ExcMessage("Object id does not correspond to an object of "
                       "the given type."));
    return at_key<T>(objects_).at(id);
};



template <class T>
bool
ObjectsContainer::
is_object (const Index &id) const
{
    const auto &objects_T = at_key<T>(objects_);
    return objects_T.find(id) != objects_T.end();
};



bool
ObjectsContainer::
is_id_present (const Index &id) const
{
    bool found = false;
    boost::fusion::for_each(objects_, [&](const auto &objects_fusion_map)
    {
        if (found)
            return;

        const auto &objects_map = objects_fusion_map.second;
        found = objects_map.find(id) != objects_map.cend();
    });

    return false;
};



IGA_NAMESPACE_CLOSE

#include <igatools/base/objects_container.inst>
