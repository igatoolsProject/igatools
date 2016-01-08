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
#include <igatools/functions/ig_function.h>


using std::shared_ptr;
using std::string;
using std::to_string;
using namespace boost::fusion;
using boost::fusion::for_each;
using std::remove_reference;
using std::dynamic_pointer_cast;
using std::const_pointer_cast;

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
insert_object(const shared_ptr<T> object,
              const bool check_present)
{
  Assert(object != nullptr, ExcNullPtr());

  auto &objects_T = at_key<T>(objects_);
  // Before inserting the object, it check_present is true it is checked
  // it the object is already present into the container, throwing an
  // exception is such case.
  AssertThrow(!check_present ||
              std::find(objects_T.cbegin(), objects_T.cend(), object) == objects_T.cend(),
              ExcNotUnique());

  objects_T.push_back(object);
};



template <class T>
void
ObjectsContainer::
insert_const_object(const shared_ptr<const T> object,
                    const bool check_present)
{
  insert_object<const T>(object,check_present);
  /*
    Assert(object != nullptr, ExcNullPtr());

    auto &objects_T = at_key<const T>(objects_);
    // Before inserting the object, it check_present is true it is checked
    // it the object is already present into the container, throwing an
    // exception is such case.
    AssertThrow(!check_present ||
                std::find(objects_T.cbegin(), objects_T.cend(), object) == objects_T.cend(),
                ExcNotUnique());

    objects_T.push_back(object);
    //*/
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
  return get_object<const T>(id);
  /*
  auto &objects_T = at_key<const T>(objects_);
  const auto obj_it = std::find_if(objects_T.cbegin(), objects_T.cend(),
  [id](const auto obj)
  {
    return obj->get_object_id() == id;
  });

  Assert(obj_it != objects_T.cend(), ExcNotFound());

  return *obj_it;
  //*/
};



template <class T>
SafeSTLSet<Index>
ObjectsContainer::
get_object_ids() const
{
  SafeSTLSet<Index> ids;
  for (const auto &obj : at_key<T>(objects_))
    ids.insert(obj->get_object_id());
  return ids;
};



template <class T>
SafeSTLSet<Index>
ObjectsContainer::
get_const_object_ids() const
{
  return get_object_ids<const T>();
  /*
    SafeSTLSet<Index> ids;
    for (const auto &obj : at_key<const T>(objects_))
      ids.insert(obj->get_object_id());
    return ids;
    //*/
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
  return is_object_present<const T>(id);
  /*
    const auto &objects_T = at_key<const T>(objects_);

    return std::find_if(objects_T.cbegin(), objects_T.cend(),
                        [id](const auto obj)
    {
      return obj->get_object_id() == id;
    }) != objects_T.cend();

    return false;
    //*/
};



void
ObjectsContainer::
print_info(LogStream &out) const
{
  auto print_info_and_id = [&out](const string &msg,const auto &objects)
  {
    out.begin_item(msg +
                   "Number of objects: " + to_string(objects.size()));
    for (const auto &obj : objects)
    {
      out.begin_item("Object Id: " + std::to_string(obj->get_object_id()));
      obj->print_info(out);
      out.end_item();
    }
    out.end_item();
  };

  GridPtrs valid_grid_ptr_types;
  for_each(valid_grid_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string msg = "Grid Dim : " + to_string(Type::dim) + ". ";
    print_info_and_id(msg,
                      at_key<Type>(objects_));

    print_info_and_id("const " + msg,
                      at_key<const Type>(objects_));
  });

  // Spline spaces
  SpSpacePtrs valid_ssp_ptr_types;
  for_each(valid_ssp_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string msg = "SplineSpace"
                       " Dim : " + to_string(Type::dim) +
                       " Range : " + to_string(Type::range) +
                       " Rank : "+ to_string(Type::rank)
                       + ". ";
    print_info_and_id(msg,
                      at_key<Type>(objects_));

    print_info_and_id("const " + msg,
                      at_key<const Type>(objects_));
  });

  // Reference space basis
  RefSpacePtrs valid_rsp_ptr_types;
  for_each(valid_rsp_ptr_types, [&](const auto &ptr_type)
  {
    using SSType = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    using Type = ReferenceSpaceBasis<SSType::dim, SSType::range, SSType::rank>;

    const string msg = "ReferenceSpaceBasis"
                       " Dim : " + to_string(Type::dim) +
                       " Range : " + to_string(Type::range) +
                       " Rank : "+ to_string(Type::rank)
                       + ". ";
    print_info_and_id(msg,
                      at_key<Type>(objects_));

    print_info_and_id("const " + msg,
                      at_key<const Type>(objects_));
  });

  // Grid functions
  GridFuncPtrs valid_gf_ptr_types;
  for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string msg = "GridFunction"
                       " Dim : " + to_string(Type::dim) +
                       " Spacedim : " + to_string(Type::range)
                       + ". ";
    print_info_and_id(msg,
                      at_key<Type>(objects_));

    print_info_and_id("const " + msg,
                      at_key<const Type>(objects_));
  });

  // Domains
  DomainPtrs valid_dm_ptr_types;
  for_each(valid_dm_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string msg = "Domain"
                       " Dim : " + to_string(Type::dim) +
                       " Codim : " + to_string(Type::space_dim - Type::dim)
                       + ". ";
    print_info_and_id(msg,
                      at_key<Type>(objects_));

    print_info_and_id("const " + msg,
                      at_key<const Type>(objects_));
  });

  // Physical space basis
  PhysSpacePtrs valid_ps_ptr_types;
  for_each(valid_ps_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string msg = "PhysicalSpaceBasis"
                       " Dim : " + to_string(Type::dim) +
                       " Range : " + to_string(Type::range) +
                       " Rank : " + to_string(Type::rank) +
                       " Codim : " + to_string(Type::codim)
                       + ". ";
    print_info_and_id(msg,
                      at_key<Type>(objects_));

    print_info_and_id("const " + msg,
                      at_key<const Type>(objects_));
  });

  // Function
  FunctionPtrs valid_fn_ptr_types;
  for_each(valid_fn_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string msg = "Function"
                       " Dim : " + to_string(Type::dim) +
                       " Codim : " + to_string(Type::codim) +
                       " Range : " + to_string(Type::range) +
                       " Rank : " + to_string(Type::rank)
                       + ". ";
    print_info_and_id(msg,
                      at_key<Type>(objects_));

    print_info_and_id("const " + msg,
                      at_key<const Type>(objects_));
  });
}



bool
ObjectsContainer::
is_empty() const
{
  auto lambda_func_not_empty = [](const auto &map_pair)
  {
    return !map_pair.second.empty();
  };

  const bool is_not_empty = boost::fusion::any(
                              objects_,
                              lambda_func_not_empty);

  return !is_not_empty;

#if 0
  bool is_container_not_void = false;

  GridPtrs valid_grid_ptr_types;
  for_each(valid_grid_ptr_types, [&](const auto &ptr_type)
  {
    if (is_container_not_void)
      return;

    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    is_container_not_void = at_key<Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;

    is_container_not_void = at_key<const Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;
  });

  if (is_container_not_void)
    return false;


  // Spline spaces
  SpSpacePtrs valid_ssp_ptr_types;
  for_each(valid_ssp_ptr_types, [&](const auto &ptr_type)
  {
    if (is_container_not_void)
      return;

    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    is_container_not_void = at_key<Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;

    is_container_not_void = at_key<const Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;
  });

  if (is_container_not_void)
    return false;


  // Reference space basis
  RefSpacePtrs valid_rsp_ptr_types;
  for_each(valid_rsp_ptr_types, [&](const auto &ptr_type)
  {
    if (is_container_not_void)
      return;

    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    is_container_not_void = at_key<Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;

    is_container_not_void = at_key<const Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;
  });

  if (is_container_not_void)
    return false;


  // Grid functions
  GridFuncPtrs valid_gf_ptr_types;
  for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
  {
    if (is_container_not_void)
      return;

    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    is_container_not_void = at_key<Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;

    is_container_not_void = at_key<const Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;
  });

  if (is_container_not_void)
    return false;


  // Domains
  DomainPtrs valid_dm_ptr_types;
  for_each(valid_dm_ptr_types, [&](const auto &ptr_type)
  {
    if (is_container_not_void)
      return;

    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    is_container_not_void = at_key<Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;

    is_container_not_void = at_key<const Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;
  });

  if (is_container_not_void)
    return false;


  // Physical space basis
  PhysSpacePtrs valid_ps_ptr_types;
  for_each(valid_ps_ptr_types, [&](const auto &ptr_type)
  {
    if (is_container_not_void)
      return;

    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    is_container_not_void = at_key<Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;

    is_container_not_void = at_key<const Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;
  });

  if (is_container_not_void)
    return false;


  // Function
  FunctionPtrs valid_fn_ptr_types;
  for_each(valid_fn_ptr_types, [&](const auto &ptr_type)
  {
    if (is_container_not_void)
      return;

    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;
    is_container_not_void = at_key<Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;

    is_container_not_void = at_key<const Type>(objects_).size() > 0;
    if (is_container_not_void)
      return;
  });

  if (is_container_not_void)
    return false;
  else
    return true;
#endif

}



void
ObjectsContainer::
fill_not_inserted_dependencies()
{
  // Adding members depending on functions (domains and physical space bases).
  FunctionPtrs valid_f_ptr_types;
  for_each(valid_f_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using DomainType = Domain<Type::dim, Type::codim>;
    using IgFunctionType = IgFunction<Type::dim, Type::codim, Type::range, Type::rank>;
    using PhysBasisType = PhysicalSpaceBasis<Type::dim, Type::range, Type::rank, Type::codim>;

    // Adding non-const objects.
    for (const auto &id : this->template get_object_ids<Type>())
    {
      // Inserting the function into the this.
      const auto obj = this->template get_object<Type>(id);
      this->template insert_object<Type> (obj);

      // Inserting the domain of the function into the this.
      const auto domain = const_pointer_cast<DomainType>(obj->get_domain());
      Assert(domain != nullptr, ExcNullPtr());
      this->template insert_object<DomainType> (domain);

      // If the function is an ig function, its physical space basis is also inserted.
      const auto const_ig_func = dynamic_pointer_cast<IgFunctionType>(obj);
      if (const_ig_func != nullptr)
      {
        const auto ig_func = const_pointer_cast<IgFunctionType>(const_ig_func);
        Assert(ig_func != nullptr, ExcNullPtr());

        const auto phys_space = const_pointer_cast<PhysBasisType>(ig_func->get_basis());
        Assert(phys_space != nullptr, ExcNullPtr());
        this->template insert_object<PhysBasisType> (phys_space);
      }
    }

    // Adding const objects.
    for (const auto &id : this->template get_const_object_ids<Type>())
    {
      // Inserting the function into the this.
      const auto obj = this->template get_const_object<Type>(id);
      this->insert_const_object<Type> (obj);

      // Inserting the domain of the function into the this.
      this->insert_const_object<DomainType> (obj->get_domain());

      // If the function is an ig function, its physical space basis is also inserted.
      const auto ig_func = dynamic_pointer_cast<const IgFunctionType>(obj);
      if (ig_func != nullptr)
        this->insert_const_object<PhysBasisType> (ig_func->get_basis());
    }
  });


  // Filling all physical space bases.
  PhysSpacePtrs valid_ps_ptr_types;
  for_each(valid_ps_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using DomainType = typename Type::PhysDomain;
    using RefBasisType = typename Type::RefBasis;

    // Adding non-const objects.
    for (const auto &id : this->template get_object_ids<Type>())
    {
      // Inserting the physical space basis into the this.
      const auto obj = this->template get_object<Type>(id);
      this->template insert_object<Type> (obj);

      // Inserting the domain of the physical space basis into the this.
      const auto domain = const_pointer_cast<DomainType>(obj->get_domain());
      Assert(domain != nullptr, ExcNullPtr());
      this->template insert_object<DomainType> (domain);

      // Inserting the reference space basis of the physical space basis into the this.
      const auto ref_space = const_pointer_cast<RefBasisType>(obj->get_reference_basis());
      Assert(ref_space != nullptr, ExcNullPtr());
      this->template insert_object<RefBasisType> (ref_space);
    }

    // Adding const objects.
    for (const auto &id : this->template get_const_object_ids<Type>())
    {
      // Inserting the physical space basis into the this.
      const auto obj = this->template get_const_object<Type>(id);
      this->insert_const_object<Type> (obj);

      // Inserting the domain of the physical space basis into the this.
      this->insert_const_object<DomainType> (obj->get_domain());

      // Inserting the reference space basis of the physical space basis into the this.
      this->insert_const_object<RefBasisType> (obj->get_reference_basis());
    }
  });


  // Filling all domains.
  DomainPtrs valid_dm_ptr_types;
  for_each(valid_dm_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridFuncType = typename Type::GridFuncType;

    // Adding non-const objects.
    for (const auto &id : this->template get_object_ids<Type>())
    {
      // Inserting the domain into the this.
      const auto obj = this->template get_object<Type>(id);
      this->template insert_object<Type> (obj);

      // Inserting the grid function of the domain into the this.
      const auto grid_func = const_pointer_cast<GridFuncType>(obj->get_grid_function());
      Assert(grid_func != nullptr, ExcNullPtr());
      this->template insert_object<GridFuncType> (grid_func);
    }

    // Adding const objects.
    for (const auto &id : this->template get_const_object_ids<Type>())
    {
      // Inserting the domain into the this.
      const auto obj = this->template get_const_object<Type>(id);
      this->insert_const_object<Type> (obj);

      // Inserting the grid function of the domain into the this.
      this->insert_const_object<GridFuncType> (obj->get_grid_function());
    }
  });


  // Filling all grid functions.
  GridFuncPtrs valid_gr_f_ptr_types;
  for_each(valid_gr_f_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridType = typename Type::GridType;
    using IgGridFuncType = IgGridFunction<Type::dim, Type::range>;
    using RefBasisType = typename IgGridFuncType::RefBasis;

    // Adding non-const objects.
    for (const auto &id : this->template get_object_ids<Type>())
    {
      // Inserting grid function into the this.
      const auto obj = this->template get_object<Type>(id);
      this->template insert_object<Type> (obj);

      // Inserting the grid of the grid function into the this.
      const auto grid = const_pointer_cast<GridType>(obj->get_grid());
      this->template insert_object<GridType> (grid);

      // If the grid function is an ig grid function, its
      // reference space basis is also inserted.
      const auto ig_g_f = dynamic_pointer_cast<IgGridFuncType>(obj);
      if (ig_g_f != nullptr)
      {
        const auto ref_space = const_pointer_cast<RefBasisType> (ig_g_f->get_basis());
        Assert(ref_space != nullptr, ExcNullPtr());
        this->template insert_object<RefBasisType> (ref_space);
      }
    }

    // Adding const objects.
    for (const auto &id : this->template get_const_object_ids<Type>())
    {
      // Inserting grid function into the this.
      const auto obj = this->template get_const_object<Type>(id);
      this->insert_const_object<Type> (obj);

      // Inserting the grid of the grid function into the this.
      this->insert_const_object<GridType> (obj->get_grid());

      // If the grid function is an ig grid function, its
      // reference space basis is also inserted.
      const auto ig_g_f = dynamic_pointer_cast<const IgGridFuncType>(obj);
      if (ig_g_f != nullptr)
        this->insert_const_object<RefBasisType> (ig_g_f->get_basis());
    }
  });


  // Filling all reference space bases.
  RefSpacePtrs valid_rs_ptr_types;
  for_each(valid_rs_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    static const int dim = Type::dim;
    static const int range = Type::range;
    static const int rank = Type::rank;
    using NURBSType = NURBS<dim, range, rank>;
    using BSplineType = typename NURBSType::BSpBasis;
    using SpSpaceType = typename BSplineType::SpSpace;
    using WeightFuncType = typename NURBSType::WeightFunction;
    using WeightFuncBasisType = typename WeightFuncType::RefBasis;
    using WeightFuncSpSpaceType = typename WeightFuncBasisType::SpSpace;
    using GridFuncType = GridFunction<WeightFuncType::dim, WeightFuncType::range>;

    // Adding non-const objects.
    for (const auto &id : this->template get_object_ids<Type>())
    {
      // Inserting reference space basis function into the this.
      const auto obj = this->template get_object<Type>(id);
      this->template insert_object<Type> (obj);

      const auto nr = dynamic_pointer_cast<NURBSType>(obj);
      const auto bs = dynamic_pointer_cast<BSplineType>(obj);

      Assert(bs != nullptr || nr != nullptr,
             ExcMessage("Invalid reference space basis type."));

      // If the reference space is a BSpline, the spline space is also inserted.
      if (bs != nullptr) // BSpline
      {
        // Inserting spline space into the this.
        const auto ss = const_pointer_cast<SpSpaceType> (bs->get_spline_space());
        this->template insert_object<SpSpaceType> (ss);
      }
      // If the reference space is a NURBS, its BSpline space, the
      // spline space, the weight grid function, the reference space of
      // the weights and spline space space of the weights are also
      // inserted.
      else // NURBS
      {
        // Inserting spline space into the this.
        const auto ss = const_pointer_cast<SpSpaceType> (nr->get_spline_space());
        this->template insert_object<SpSpaceType> (ss);

        // Inserting the BSpline space of the NURBS into the this.
        const auto nr_bs = const_pointer_cast<BSplineType> (nr->get_bspline_basis());
        this->template insert_object<Type> (nr_bs);

        // Adding the weight related quantities: grid function,
        // reference space basis and its spline space.

        // Inserting the grid function of the weights.
        const auto wf = const_pointer_cast<WeightFuncType>(nr->get_weight_func());
        Assert(wf != nullptr, ExcNullPtr());
        this->template insert_object<GridFuncType>(wf);

        // Inserting the reference space basis of the weights.
        const auto rs = const_pointer_cast<WeightFuncBasisType>(wf->get_basis());
        Assert(rs != nullptr, ExcNullPtr());
        this->template insert_object<WeightFuncBasisType>(rs);

        // Inserting the spline space basis of the weights.
        const auto rs_ss = const_pointer_cast<WeightFuncSpSpaceType>(rs->get_spline_space());
        this->template insert_object<WeightFuncSpSpaceType>(rs_ss);
      }
    }

    // Adding const objects.
    for (const auto &id : this->template get_const_object_ids<Type>())
    {
      // Inserting reference space basis function into the this.
      const auto obj = this->template get_const_object<Type>(id);
      this->insert_const_object<Type> (obj);

      const auto nr = dynamic_pointer_cast<const NURBSType>(obj);
      const auto bs = dynamic_pointer_cast<const BSplineType>(obj);

      Assert(bs != nullptr || nr != nullptr,
             ExcMessage("Invalid reference space basis type."));

      // If the reference space is a BSpline, the spline space is also inserted.
      if (bs != nullptr) // BSpline
      {
        // Inserting spline space into the this.
        this->insert_const_object<SpSpaceType> (bs->get_spline_space());
      }
      // If the reference space is a NURBS, its BSpline space, the
      // spline space, the weight grid function, the reference space of
      // the weights and spline space space of the weights are also
      // inserted.
      else // NURBS
      {
        // Inserting spline space into the this.
        this->insert_const_object<SpSpaceType> (nr->get_spline_space());

        // Inserting the BSpline space of the NURBS into the this.
        this->insert_const_object<Type> (nr->get_bspline_basis());

        // Adding the weight related quantities: grid function,
        // reference space basis and its spline space.

        const auto wf = nr->get_weight_func();
        const auto rs = wf->get_basis();
        const auto gf = dynamic_pointer_cast<const GridFuncType> (wf);
        Assert(gf != nullptr, ExcNullPtr());

        // Inserting the grid function of the weights.
        this->insert_const_object<GridFuncType>(gf);

        // Inserting the reference space basis of the weights.
        this->insert_const_object<WeightFuncBasisType>(rs);

        // Inserting the spline space basis of the weights.
        this->insert_const_object<WeightFuncSpSpaceType>(rs->get_spline_space());
      }
    }
  });


  // Filling all spline spaces.
  SpSpacePtrs valid_ss_ptr_types;
  for_each(valid_ss_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridType = Grid<Type::dim>;

    // Adding non-const objects.
    for (const auto &id : this->template get_object_ids<Type>())
    {
      // Inserting spline space into the this.
      const auto obj = this->template get_object<Type>(id);
      this->template insert_object<Type> (obj);

      // Inserting the grid of the spline space into the this.
      const auto grid = const_pointer_cast<GridType>(obj->get_grid());
      this->template insert_object<GridType> (grid);
    }

    // Adding const objects.
    for (const auto &id : this->template get_const_object_ids<Type>())
    {
      // Inserting spline space into the this.
      const auto obj = this->template get_const_object<Type>(id);
      this->insert_const_object<Type> (obj);

      // Inserting the grid of the spline space into the this.
      this->insert_const_object<GridType> (obj->get_grid());
    }
  });

  // Filling all grids.
  GridPtrs valid_gr_ptr_types;
  for_each(valid_gr_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    // Adding non-const objects.
    for (const auto &id : this->template get_object_ids<Type>())
    {
      // Inserting grid into the this.
      const auto obj = this->template get_object<Type>(id);
      this->template insert_object<Type> (obj);
    }

    // Adding const objects.
    for (const auto &id : this->template get_const_object_ids<Type>())
    {
      // Inserting grid into the this.
      const auto obj = this->template get_const_object<Type>(id);
      this->insert_const_object<Type> (obj);
    }
  });
}

#ifdef SERIALIZATION

template<class Archive>
void
ObjectsContainer::
serialize(Archive &ar)
{
  // Serializing grids.
  GridPtrs valid_grid_ptr_types;
  for_each(valid_grid_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string name = "grid_"  + to_string(Type::dim);

    ar &make_nvp(name, at_key<Type>(objects_));

    auto &const_objects = at_key<const Type>(objects_);

    SafeSTLVector<shared_ptr<Type>> tmp_objects;
    for (auto &obj : const_objects)
      tmp_objects.push_back(const_pointer_cast<Type>(obj));

    ar &make_nvp("const_" + name, tmp_objects);

    if (const_objects.empty())
      for (const auto &obj : tmp_objects)
        const_objects.push_back(obj);
  });

  // Serializing spline spaces
  SpSpacePtrs valid_ssp_ptr_types;
  for_each(valid_ssp_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string name = "spline_space_"
                        + to_string(Type::dim) + "_"
                        + to_string(Type::range) + "_"
                        + to_string(Type::rank);

    ar &make_nvp(name, at_key<Type>(objects_));
  });

  // Serializing reference space basis
  RefSpacePtrs valid_rsp_ptr_types;
  for_each(valid_rsp_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string name = "reference_space_basis_"
                        + to_string(Type::dim) + "_"
                        + to_string(Type::range) + "_"
                        + to_string(Type::rank);

    ar &make_nvp(name, at_key<Type>(objects_));

    auto &const_objects = at_key<const Type>(objects_);

    SafeSTLVector<shared_ptr<Type>> tmp_objects;
    for (auto &obj : const_objects)
      tmp_objects.push_back(const_pointer_cast<Type>(obj));

    ar &make_nvp("const_" + name, tmp_objects);

    if (const_objects.empty())
      for (const auto &obj : tmp_objects)
        const_objects.push_back(obj);
  });

  // Grid functions
  GridFuncPtrs valid_gf_ptr_types;
  for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string name = "grid_funcion_"
                        + to_string(Type::dim) + "_"
                        + to_string(Type::range);

    ar &make_nvp(name, at_key<Type>(objects_));

    auto &const_objects = at_key<const Type>(objects_);

    SafeSTLVector<shared_ptr<Type>> tmp_objects;
    for (auto &obj : const_objects)
      tmp_objects.push_back(const_pointer_cast<Type>(obj));

    ar &make_nvp("const_" + name, tmp_objects);

    if (const_objects.empty())
      for (const auto &obj : tmp_objects)
        const_objects.push_back(obj);
  });

  // Domains
  DomainPtrs valid_dm_ptr_types;
  for_each(valid_dm_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string name = "domain_"
                        + to_string(Type::dim) + "_"
                        + to_string(Type::space_dim);

    ar &make_nvp(name, at_key<Type>(objects_));

    auto &const_objects = at_key<const Type>(objects_);

    SafeSTLVector<shared_ptr<Type>> tmp_objects;
    for (auto &obj : const_objects)
      tmp_objects.push_back(const_pointer_cast<Type>(obj));

    ar &make_nvp("const_" + name, tmp_objects);

    if (const_objects.empty())
      for (const auto &obj : tmp_objects)
        const_objects.push_back(obj);
  });

  // Physical space basis
  PhysSpacePtrs valid_ps_ptr_types;
  for_each(valid_ps_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string name = "physical_space_basis_"
                        + to_string(Type::dim) + "_"
                        + to_string(Type::range) + "_"
                        + to_string(Type::rank) + "_"
                        + to_string(Type::codim);

    ar &make_nvp(name, at_key<Type>(objects_));

    auto &const_objects = at_key<const Type>(objects_);

    SafeSTLVector<shared_ptr<Type>> tmp_objects;
    for (auto &obj : const_objects)
      tmp_objects.push_back(const_pointer_cast<Type>(obj));

    ar &make_nvp("const_" + name, tmp_objects);

    if (const_objects.empty())
      for (const auto &obj : tmp_objects)
        const_objects.push_back(obj);
  });

  // Function
  FunctionPtrs valid_fn_ptr_types;
  for_each(valid_fn_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename std::remove_reference<decltype(ptr_type)>::type::element_type;

    const string name = "function_"
                        + to_string(Type::dim) + "_"
                        + to_string(Type::codim) + "_"
                        + to_string(Type::range) + "_"
                        + to_string(Type::rank);

    ar &make_nvp(name, at_key<Type>(objects_));

    auto &const_objects = at_key<const Type>(objects_);

    SafeSTLVector<shared_ptr<Type>> tmp_objects;
    for (auto &obj : const_objects)
      tmp_objects.push_back(const_pointer_cast<Type>(obj));

    ar &make_nvp("const_" + name, tmp_objects);

    if (const_objects.empty())
      for (const auto &obj : tmp_objects)
        const_objects.push_back(obj);
  });



}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/base/objects_container.inst>


#ifdef SERIALIZATION
template void iga::ObjectsContainer::serialize(OArchive &);
template void iga::ObjectsContainer::serialize(IArchive &);
#endif // SERIALIZATION

#endif // XML_IO
