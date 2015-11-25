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
#include <igatools/functions/function.h>

#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/container/map/map_fwd.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <boost/fusion/include/at_key.hpp>


using std::shared_ptr;
using std::to_string;
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



void
ObjectsContainer::
print_info (LogStream &out) const
{
    using GPtrs = typename InstantiatedTypes::ValidGridPtrs;
    using SSPtrs = typename InstantiatedTypes::ValidSplineSpacePtrs;
    using GFPtrs = typename InstantiatedTypes::ValidGridFunctionPtrs;
    using DPtrs = typename InstantiatedTypes::ValidDomainPtrs;
    using PSPtrs = typename InstantiatedTypes::ValidPhysSpacePtrs;
    using FPtrs = typename InstantiatedTypes::ValidFunctionPtrs;


    // Grids
    GPtrs valid_grid_ptr_types;
    boost::fusion::for_each(valid_grid_ptr_types, [&](const auto &ptr_type)
    {
        using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
        const auto objects = at_key<Type>(objects_);

        out.begin_item("Grid"
                " Dim : " + to_string(Type::dim) +
                ". Number of objects: " + to_string(objects.size()));
        for (const auto &object : objects)
        {
            out.begin_item("Object Id : " + to_string(object.first));
            object.second->print_info(out);
            out.end_item();
        }
        out.end_item();
    });

    // Spline spaces
    SSPtrs valid_ssp_ptr_types;
    boost::fusion::for_each(valid_ssp_ptr_types, [&](const auto &ptr_type)
    {
        using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
        const auto objects = at_key<Type>(objects_);

        out.begin_item("SplineSpace"
                " Dim : " + to_string(Type::dim) +
                " Range : " + to_string(Type::range) +
                " Rank : "+ to_string(Type::rank) +
                ". Number of objects: " + to_string(objects.size())
        );
        for (const auto &object : objects)
        {
            out.begin_item("Object Id : " + to_string(object.first));
            object.second->print_info(out);
            out.end_item();
        }
        out.end_item();
    });

    // Reference space basis
    boost::fusion::for_each(valid_ssp_ptr_types, [&](const auto &ptr_type)
    {
        using SSType = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
        using Type = ReferenceSpaceBasis<SSType::dim, SSType::range, SSType::rank>;
        const auto objects = at_key<Type>(objects_);

        out.begin_item("ReferenceSpaceBasis"
                " Dim : " + to_string(Type::dim) +
                " Range : " + to_string(Type::range) +
                " Rank : "+ to_string(Type::rank) +
                ". Number of objects: " + to_string(objects.size())
        );
        for (const auto &object : objects)
        {
            out.begin_item("Object Id : " + to_string(object.first));
            object.second->print_info(out);
            out.end_item();
        }
        out.end_item();
    });

    // Grid functions
    GFPtrs valid_gf_ptr_types;
    boost::fusion::for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
    {
        using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
        const auto objects = at_key<Type>(objects_);

        out.begin_item("GridFunction"
                " Dim : " + to_string(Type::dim) +
                " Spacedim : " + to_string(Type::space_dim) +
                ". Number of objects: " + to_string(objects.size()));
        for (const auto &object : objects)
        {
            out.begin_item("Object Id : " + to_string(object.first));
            object.second->print_info(out);
            out.end_item();
        }
        out.end_item();
    });

    // Domains
    DPtrs valid_dm_ptr_types;
    boost::fusion::for_each(valid_dm_ptr_types, [&](const auto &ptr_type)
    {
        using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
        const auto objects = at_key<Type>(objects_);

        out.begin_item("Domain"
                " Dim : " + to_string(Type::dim) +
                " Codim : " + to_string(Type::space_dim - Type::dim) +
                ". Number of objects: " + to_string(objects.size()));
        for (const auto &object : objects)
        {
            out.begin_item("Object Id : " + to_string(object.first));
            object.second->print_info(out);
            out.end_item();
        }
        out.end_item();
    });

    // Physical space basis
    PSPtrs valid_ps_ptr_types;
    boost::fusion::for_each(valid_ps_ptr_types, [&](const auto &ptr_type)
    {
        using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
        const auto objects = at_key<Type>(objects_);

        out.begin_item("PhysicalSpaceBasis"
                " Dim : " + to_string(Type::dim) +
                " Range : " + to_string(Type::range) +
                " Rank : " + to_string(Type::rank) +
                " Codim : " + to_string(Type::codim) +
                ". Number of objects: " + to_string(objects.size()));
        for (const auto &object : objects)
        {
            out.begin_item("Object Id : " + to_string(object.first));
            object.second->print_info(out);
            out.end_item();
        }
        out.end_item();
    });

    // Function
    FPtrs valid_fn_ptr_types;
    boost::fusion::for_each(valid_fn_ptr_types, [&](const auto &ptr_type)
    {
        using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
        const auto objects = at_key<Type>(objects_);

        out.begin_item("Function"
                " Dim : " + to_string(Type::dim) +
                " Codim : " + to_string(Type::codim) +
                " Range : " + to_string(Type::range) +
                " Rank : " + to_string(Type::rank) +
                ". Number of objects: " + to_string(objects.size()));
        for (const auto &object : objects)
        {
            out.begin_item("Object Id : " + to_string(object.first));
            object.second->print_info(out);
            out.end_item();
        }
        out.end_item();
    });
}



IGA_NAMESPACE_CLOSE

#include <igatools/base/objects_container.inst>
