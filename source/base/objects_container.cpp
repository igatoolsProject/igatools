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

#ifdef XML_IO

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/nurbs.h>
#include <igatools/functions/function.h>


using std::shared_ptr;
using std::to_string;
using namespace boost::fusion;
using boost::fusion::for_each;

IGA_NAMESPACE_OPEN

auto
ObjectsContainer::
create() -> shared_ptr<self_t>
{
  return shared_ptr<self_t> (new self_t());
};



template <class T>
void
ObjectsContainer::
insert_object(const shared_ptr<T> object)
{
  Assert(object != nullptr, ExcNullPtr());

  auto &objects_T = at_key<T>(objects_);
  Assert(std::find(objects_T.cbegin(), objects_T.cend(), object) == objects_T.cend(),
         ExcNotUnique());

  objects_T.push_back(object);
};



template <class T>
void
ObjectsContainer::
insert_const_object(const shared_ptr<const T> object)
{
  Assert(object != nullptr, ExcNullPtr());

  auto &objects_T = at_key<const T>(objects_);
  Assert(std::find(objects_T.cbegin(), objects_T.cend(), object) == objects_T.cend(),
         ExcNotUnique());

  objects_T.push_back(object);
};



template <class T>
auto
ObjectsContainer::
get_object(const Index &id) const -> shared_ptr<T>
{
  auto &objects_T = at_key<T>(objects_);
  const auto obj_it = std::find_if(objects_T.cbegin(), objects_T.cend(),
  [id](const auto obj)
  {
    return obj->get_object_id() == id;
  });

  Assert(obj_it != objects_T.cend(), ExcNotFound());

  return *obj_it;
};



template <class T>
auto
ObjectsContainer::
get_const_object(const Index &id) const -> shared_ptr<const T>
{
  auto &objects_T = at_key<const T>(objects_);
  const auto obj_it = std::find_if(objects_T.cbegin(), objects_T.cend(),
  [id](const auto obj)
  {
    return obj->get_object_id() == id;
  });

  Assert(obj_it != objects_T.cend(), ExcNotFound());

  return *obj_it;
};



template <class T>
bool
ObjectsContainer::
is_object_present(const Index &id) const
{
  const auto &objects_T = at_key<T>(objects_);

  return std::find_if(objects_T.cbegin(), objects_T.cend(),
                      [id](const auto obj)
  {
    return obj->get_object_id() == id;
  }) != objects_T.cend();

  return false;
};



template <class T>
bool
ObjectsContainer::
is_const_object_present(const Index &id) const
{
  const auto &objects_T = at_key<const T>(objects_);

  return std::find_if(objects_T.cbegin(), objects_T.cend(),
                      [id](const auto obj)
  {
    return obj->get_object_id() == id;
  }) != objects_T.cend();

  return false;
};



void
ObjectsContainer::
print_info(LogStream &out) const
{
  GridPtrs valid_grid_ptr_types;
  for_each(valid_grid_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    const auto objects = at_key<Type>(objects_);

    out.begin_item("Grid"
                   " Dim : " + to_string(Type::dim) +
                   ". Number of objects: " + to_string(objects.size()));
    for (const auto &object : objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();

    const auto const_objects = at_key<const Type>(objects_);

    out.begin_item("const Grid"
                   " Dim : " + to_string(Type::dim) +
                   ". Number of objects: " + to_string(objects.size()));
    for (const auto &object : const_objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();
  });

  // Spline spaces
  SpSpacePtrs valid_rsp_ptr_types;
  for_each(valid_rsp_ptr_types, [&](const auto &ptr_type)
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
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();

    const auto const_objects = at_key<const Type>(objects_);

    out.begin_item("const SplineSpace"
                   " Dim : " + to_string(Type::dim) +
                   " Range : " + to_string(Type::range) +
                   " Rank : "+ to_string(Type::rank) +
                   ". Number of objects: " + to_string(objects.size())
                  );
    for (const auto &object : const_objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();
  });

  // Reference space basis
  for_each(valid_rsp_ptr_types, [&](const auto &ptr_type)
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
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();

    const auto const_objects = at_key<const Type>(objects_);

    out.begin_item("const ReferenceSpaceBasis"
                   " Dim : " + to_string(Type::dim) +
                   " Range : " + to_string(Type::range) +
                   " Rank : "+ to_string(Type::rank) +
                   ". Number of objects: " + to_string(objects.size())
                  );
    for (const auto &object : const_objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();
  });

  // Grid functions
  GridFuncPtrs valid_gf_ptr_types;
  for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    const auto objects = at_key<Type>(objects_);

    out.begin_item("GridFunction"
                   " Dim : " + to_string(Type::dim) +
                   " Spacedim : " + to_string(Type::space_dim) +
                   ". Number of objects: " + to_string(objects.size()));
    for (const auto &object : objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();

    const auto const_objects = at_key<const Type>(objects_);

    out.begin_item("const GridFunction"
                   " Dim : " + to_string(Type::dim) +
                   " Spacedim : " + to_string(Type::space_dim) +
                   ". Number of objects: " + to_string(objects.size()));
    for (const auto &object : const_objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();
  });

  // Domains
  DomainPtrs valid_dm_ptr_types;
  for_each(valid_dm_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    const auto objects = at_key<Type>(objects_);

    out.begin_item("Domain"
                   " Dim : " + to_string(Type::dim) +
                   " Codim : " + to_string(Type::space_dim - Type::dim) +
                   ". Number of objects: " + to_string(objects.size()));
    for (const auto &object : objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();

    const auto const_objects = at_key<const Type>(objects_);

    out.begin_item("const Domain"
                   " Dim : " + to_string(Type::dim) +
                   " Codim : " + to_string(Type::space_dim - Type::dim) +
                   ". Number of objects: " + to_string(objects.size()));
    for (const auto &object : const_objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();
  });

  // Physical space basis
  PhysSpacePtrs valid_ps_ptr_types;
  for_each(valid_ps_ptr_types, [&](const auto &ptr_type)
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
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();

    const auto const_objects = at_key<const Type>(objects_);

    out.begin_item("const PhysicalSpaceBasis"
                   " Dim : " + to_string(Type::dim) +
                   " Range : " + to_string(Type::range) +
                   " Rank : " + to_string(Type::rank) +
                   " Codim : " + to_string(Type::codim) +
                   ". Number of objects: " + to_string(objects.size()));
    for (const auto &object : const_objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();
  });

  // Function
  FunctionPtrs valid_fn_ptr_types;
  for_each(valid_fn_ptr_types, [&](const auto &ptr_type)
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
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();

    const auto const_objects = at_key<const Type>(objects_);

    out.begin_item("const Function"
                   " Dim : " + to_string(Type::dim) +
                   " Codim : " + to_string(Type::codim) +
                   " Range : " + to_string(Type::range) +
                   " Rank : " + to_string(Type::rank) +
                   ". Number of objects: " + to_string(objects.size()));
    for (const auto &object : const_objects)
    {
      out.begin_item("Object Id: " + std::to_string(object->get_object_id()));
      object->print_info(out);
      out.end_item();
    }
    out.end_item();
  });
}


IGA_NAMESPACE_CLOSE

#include <igatools/base/objects_container.inst>

#endif // XML_IO
