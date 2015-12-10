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

#include <igatools/io/objects_container_writer.h>

#ifdef XML_IO

#include <igatools/base/objects_container.h>
#include <igatools/utils/safe_stl_set.h>
//#include <igatools/io/objects_container_parser-XML_schema.h>
//
//#include <igatools/io/xml_file_parser.h>
#include <igatools/io/xml_element.h>

#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_function_lib.h>
#include <igatools/basis_functions/spline_space.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/nurbs.h>
#include <igatools/functions/function.h>
#include <igatools/functions/ig_function.h>
#include <igatools/functions/function_lib.h>

#include <boost/fusion/algorithm/iteration/for_each.hpp>

using std::string;
//using std::to_string;
using std::shared_ptr;
//using std::set;
using std::remove_reference;
using std::dynamic_pointer_cast;
using std::const_pointer_cast;
//using std::copy;
//using std::inserter;

IGA_NAMESPACE_OPEN

void
ObjectsContainerWriter::
write(const string &file_path,
      const shared_ptr<ObjectsContainer> container)
{
    AssertThrow (false, ExcNotImplemented());
    const auto full_container = Self_::build_full_container(container);

    xercesc::DOMElement *dom_elem = nullptr;
    const auto xml_elem = XMLElement::create(dom_elem);

    Self_::write_grids(full_container, *xml_elem);
    Self_::write_spline_spaces(full_container, *xml_elem);
    Self_::write_reference_space_bases(full_container, *xml_elem);
    Self_::write_grid_functions(full_container, *xml_elem);
    Self_::write_domains(full_container, *xml_elem);
    Self_::write_physical_space_bases(full_container, *xml_elem);
    Self_::write_functions(full_container, *xml_elem);

    // TODO: write xml_elem to file.
}



void
ObjectsContainerWriter::
write_grids (const shared_ptr<ObjectsContainer> container,
             XMLElement &xml_elem)
{
  using GridPtrs = typename ObjectsContainer::GridPtrs;

  GridPtrs valid_grid_ptr_types;
  for_each(valid_grid_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    const auto non_const_ids = container->template get_object_ids<Type>();

    for (const auto &id : non_const_ids)
    {
        const auto grid = container->template get_object<Type>(id);
        Self_::write_grid<Type>(grid, xml_elem);
    }


    // Adding non-const objects.
    const auto const_ids = container->template get_const_object_ids<Type>();
    for (const auto &id : const_ids)
    {
        const auto grid = container->template get_const_object<Type>(id);
        Self_::write_grid<const Type>(grid, xml_elem);
    }
  });
}



void
ObjectsContainerWriter::
write_spline_spaces (const shared_ptr<ObjectsContainer> container,
             XMLElement &xml_elem)
{
  using SpSpacePtrs = typename ObjectsContainer::SpSpacePtrs;

  SpSpacePtrs valid_spline_space_ptr_types;
  for_each(valid_spline_space_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    const auto non_const_ids = container->template get_object_ids<Type>();

    for (const auto &id : non_const_ids)
    {
        const auto sp_space = container->template get_object<Type>(id);
        Self_::write_spline_space<Type>(sp_space, xml_elem);
    }


    const auto const_ids = container->template get_const_object_ids<Type>();
    for (const auto &id : const_ids)
    {
        const auto sp_space = container->template get_const_object<Type>(id);
        Self_::write_spline_space<const Type>(sp_space, xml_elem);
    }
  });
}



void
ObjectsContainerWriter::
write_reference_space_bases (const shared_ptr<ObjectsContainer> container,
             XMLElement &xml_elem)
{
  using RefSpacePtrs = typename ObjectsContainer::RefSpacePtrs;

  RefSpacePtrs valid_ref_space_ptr_types;
  for_each(valid_ref_space_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    static const int dim = Type::dim;
    static const int range = Type::range;
    static const int rank = Type::rank;

    using BSplineType = BSpline<dim, range, rank>;
    using NURBSType = NURBS<dim, range, rank>;

    const auto non_const_ids = container->template get_object_ids<Type>();

    for (const auto &id : non_const_ids)
    {
        const auto ref_space = container->template get_object<Type>(id);

        const auto bs_space = dynamic_pointer_cast<BSplineType>(ref_space);
        const auto nr_space = dynamic_pointer_cast<NURBSType>(ref_space);

        Assert (bs_space == nullptr && nr_space == nullptr,
                ExcMessage("Invalid reference space type."));

        if (bs_space != nullptr)
            Self_::write_bspline<BSplineType>(bs_space, xml_elem);
        else
            Self_::write_bspline<NURBSType>(nr_space, xml_elem);
    }


    const auto const_ids = container->template get_const_object_ids<Type>();
    for (const auto &id : const_ids)
    {
        const auto ref_space = container->template get_const_object<Type>(id);

        const auto bs_space = dynamic_pointer_cast<const BSplineType>(ref_space);
        const auto nr_space = dynamic_pointer_cast<const NURBSType>(ref_space);

        Assert (bs_space == nullptr && nr_space == nullptr,
                ExcMessage("Invalid reference space type."));

        if (bs_space != nullptr)
            Self_::write_bspline<const BSplineType>(bs_space, xml_elem);
        else
            Self_::write_bspline<const NURBSType>(nr_space, xml_elem);
    }
  });
}



void
ObjectsContainerWriter::
write_grid_functions (const shared_ptr<ObjectsContainer> container,
                      XMLElement &xml_elem)
{
  using GridFuncPtrs = typename ObjectsContainer::GridFuncPtrs;

  GridFuncPtrs valid_gf_ptr_types;
  for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    static const int dim = Type::dim;
    static const int space_dim = Type::space_dim;

    using IdGridFunc = grid_functions::IdentityGridFunction<dim>;
    using LinearGridFunc = grid_functions::LinearGridFunction<dim, space_dim>;
    using ConstantGridFunc = grid_functions::ConstantGridFunction<dim, space_dim>;
    using IgGridFunc = IgGridFunction<dim, space_dim>;

    const auto non_const_ids = container->template get_object_ids<Type>();

    for (const auto &id : non_const_ids)
    {
        const auto grid_func = container->template get_object<Type>(id);

        const auto id_f = dynamic_pointer_cast<IdGridFunc>(grid_func);
        const auto li_f = dynamic_pointer_cast<LinearGridFunc>(grid_func);
        const auto ct_f = dynamic_pointer_cast<ConstantGridFunc>(grid_func);
        const auto ig_f = dynamic_pointer_cast<IgGridFunc>(grid_func);

        Assert (id_f == nullptr && li_f == nullptr &&
                ct_f == nullptr && ig_f == nullptr,
                ExcMessage("Invalid grid function type."));

        if (id_f != nullptr)
            Self_::write_identity_grid_function<IdGridFunc>(id_f, xml_elem);
        else if (li_f != nullptr)
            Self_::write_linear_grid_function<LinearGridFunc>(li_f, xml_elem);
        else if (ct_f != nullptr)
            Self_::write_constant_grid_function<ConstantGridFunc>(ct_f, xml_elem);
        else
            Self_::write_ig_grid_function<IgGridFunc>(ig_f, xml_elem);
    }


    const auto const_ids = container->template get_const_object_ids<Type>();
    for (const auto &id : const_ids)
    {
        const auto grid_func = container->template get_const_object<Type>(id);

        const auto id_f = dynamic_pointer_cast<const IdGridFunc>(grid_func);
        const auto li_f = dynamic_pointer_cast<const LinearGridFunc>(grid_func);
        const auto ct_f = dynamic_pointer_cast<const ConstantGridFunc>(grid_func);
        const auto ig_f = dynamic_pointer_cast<const IgGridFunc>(grid_func);

        Assert (id_f == nullptr && li_f == nullptr &&
                ct_f == nullptr && ig_f == nullptr,
                ExcMessage("Invalid grid function type."));

        if (id_f != nullptr)
            Self_::write_identity_grid_function<const IdGridFunc>(id_f, xml_elem);
        else if (li_f != nullptr)
            Self_::write_linear_grid_function<const LinearGridFunc>(li_f, xml_elem);
        else if (ct_f != nullptr)
            Self_::write_constant_grid_function<const ConstantGridFunc>(ct_f, xml_elem);
        else
            Self_::write_ig_grid_function<const IgGridFunc>(ig_f, xml_elem);
    }
  });
}



void
ObjectsContainerWriter::
write_domains (const shared_ptr<ObjectsContainer> container,
               XMLElement &xml_elem)
{
  using DomainPtrs = typename ObjectsContainer::DomainPtrs;

  DomainPtrs valid_dm_ptr_types;
  for_each(valid_dm_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    const auto non_const_ids = container->template get_object_ids<Type>();

    for (const auto &id : non_const_ids)
    {
        const auto domain = container->template get_object<Type>(id);
        Self_::write_domain<Type>(domain, xml_elem);
    }


    const auto const_ids = container->template get_const_object_ids<Type>();
    for (const auto &id : const_ids)
    {
        const auto domain = container->template get_const_object<Type>(id);
        Self_::write_domain<const Type>(domain, xml_elem);
    }
  });
}



void
ObjectsContainerWriter::
write_physical_space_bases (const shared_ptr<ObjectsContainer> container,
                            XMLElement &xml_elem)
{
  using PhysSpacePtrs = typename ObjectsContainer::PhysSpacePtrs;

  PhysSpacePtrs valid_ps_ptr_types;
  for_each(valid_ps_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    const auto non_const_ids = container->template get_object_ids<Type>();

    for (const auto &id : non_const_ids)
    {
        const auto ps_space = container->template get_object<Type>(id);
        Self_::write_phys_space_basis<Type>(ps_space, xml_elem);
    }


    const auto const_ids = container->template get_const_object_ids<Type>();
    for (const auto &id : const_ids)
    {
        const auto ps_space = container->template get_const_object<Type>(id);
        Self_::write_phys_space_basis<const Type>(ps_space, xml_elem);
    }
  });
}



void
ObjectsContainerWriter::
write_functions (const shared_ptr<ObjectsContainer> container,
                 XMLElement &xml_elem)
{
  using FuncPtrs = typename ObjectsContainer::FunctionPtrs;

  FuncPtrs valid_f_ptr_types;
  for_each(valid_f_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    static const int dim = Type::dim;
    static const int codim = Type::codim;
    static const int range = Type::range;
    static const int rank = Type::rank;

    using LinearFunc = functions::LinearFunction<dim, codim, range>;
    using ConstantFunc = functions::ConstantFunction<dim, codim, range, rank>;
    using IgFunc = IgFunction<dim, codim, range, rank>;

    const auto non_const_ids = container->template get_object_ids<Type>();

    for (const auto &id : non_const_ids)
    {
        const auto func = container->template get_object<Type>(id);

        const auto li_f = dynamic_pointer_cast<LinearFunc>(func);
        const auto ct_f = dynamic_pointer_cast<ConstantFunc>(func);
        const auto ig_f = dynamic_pointer_cast<IgFunc>(func);

        Assert (li_f == nullptr && ct_f == nullptr && ig_f == nullptr,
                ExcMessage("Invalid function type."));

        if (li_f != nullptr)
            Self_::write_linear_function<LinearFunc>(li_f, xml_elem);
        else if (ct_f != nullptr)
            Self_::write_constant_function<ConstantFunc>(ct_f, xml_elem);
        else
            Self_::write_ig_function<IgFunc>(ig_f, xml_elem);
    }


    const auto const_ids = container->template get_const_object_ids<Type>();
    for (const auto &id : const_ids)
    {
        const auto func = container->template get_const_object<Type>(id);

        const auto li_f = dynamic_pointer_cast<const LinearFunc>(func);
        const auto ct_f = dynamic_pointer_cast<const ConstantFunc>(func);
        const auto ig_f = dynamic_pointer_cast<const IgFunc>(func);

        Assert (li_f == nullptr && ct_f == nullptr && ig_f == nullptr,
                ExcMessage("Invalid function type."));

        if (li_f != nullptr)
            Self_::write_linear_function<const LinearFunc>(li_f, xml_elem);
        else if (ct_f != nullptr)
            Self_::write_constant_function<const ConstantFunc>(ct_f, xml_elem);
        else
            Self_::write_ig_function<const IgFunc>(ig_f, xml_elem);
    }
  });
}



template <class Grid>
void
ObjectsContainerWriter::
write_grid (const shared_ptr<Grid> grid,
            XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class SpSpace>
void
ObjectsContainerWriter::
write_spline_space (const shared_ptr<SpSpace> spline_space,
                    XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class BSpline>
void
ObjectsContainerWriter::
write_bspline (const shared_ptr<BSpline> bspline,
               XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class NURBS>
void
ObjectsContainerWriter::
write_nurbs (const shared_ptr<NURBS> nurbs,
             XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class IdGridFunc>
void
ObjectsContainerWriter::
write_identity_grid_function (const shared_ptr<IdGridFunc> id_func,
                              XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class ConstGridFunc>
void
ObjectsContainerWriter::
write_constant_grid_function (const shared_ptr<ConstGridFunc> const_func,
                              XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class LinearGridFunc>
void
ObjectsContainerWriter::
write_linear_grid_function (const shared_ptr<LinearGridFunc> linear_func,
                            XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class IgGridFunc>
void
ObjectsContainerWriter::
write_ig_grid_function (const shared_ptr<IgGridFunc> ig_func,
                        XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class Domain>
void
ObjectsContainerWriter::
write_domain (const shared_ptr<Domain> domain,
              XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class PhysSpaceBasis>
void
ObjectsContainerWriter::
write_phys_space_basis (const shared_ptr<PhysSpaceBasis> phys_space,
                        XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class IgFunction>
void
ObjectsContainerWriter::
write_ig_function (const shared_ptr<IgFunction> ig_function,
                   XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class ConstantFunction>
void
ObjectsContainerWriter::
write_constant_function (const shared_ptr<ConstantFunction> const_function,
                         XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class LinearFunction>
void
ObjectsContainerWriter::
write_linear_function (const shared_ptr<LinearFunction> linear_function,
                         XMLElement &xml_elem)
{
    AssertThrow (false, ExcNotImplemented());
}



shared_ptr<ObjectsContainer>
ObjectsContainerWriter::
build_full_container (const shared_ptr<ObjectsContainer> container)
{
  const auto full_container = ObjectsContainer::create ();

  using FuncPtrs = typename ObjectsContainer::FunctionPtrs;
  using PhysSpacePtrs = typename ObjectsContainer::PhysSpacePtrs;
  using DomainPtrs = typename ObjectsContainer::DomainPtrs;
  using GridFuncPtrs = typename ObjectsContainer::GridFuncPtrs;
  using RefSpacePtrs = typename ObjectsContainer::RefSpacePtrs;
  using SpSpacePtrs = typename ObjectsContainer::SpSpacePtrs;
  using GridPtrs = typename ObjectsContainer::GridPtrs;

  // Adding members depending on functions (domains and physical space bases).
  FuncPtrs valid_f_ptr_types;
  for_each(valid_f_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using DomainType = Domain<Type::dim, Type::codim>;
    using IgFunctionType = IgFunction<Type::dim, Type::codim, Type::range, Type::rank>;
    using PhysBasisType = PhysicalSpaceBasis<Type::dim, Type::range, Type::rank, Type::codim>;

    // Adding const objects.
    for (const auto &id : container->template get_object_ids<Type>())
    {
        // Inserting the function into the container.
        const auto obj = container->template get_object<Type>(id);
        full_container->insert_object<Type> (obj);

        // Inserting the domain of the function into the container.
        const auto domain = const_pointer_cast<DomainType>(obj->get_domain());
        Assert (domain != nullptr, ExcNullPtr());
        full_container->insert_object<DomainType> (domain);

        // If the function is an ig function, its physical space basis is also inserted.
        const auto const_ig_func = dynamic_pointer_cast<IgFunctionType>(obj);
        if (const_ig_func != nullptr)
        {
            const auto ig_func = const_pointer_cast<IgFunctionType>(const_ig_func);
            Assert (ig_func != nullptr, ExcNullPtr());

            const auto phys_space = const_pointer_cast<PhysBasisType>(ig_func->get_basis());
            Assert (phys_space != nullptr, ExcNullPtr());
            full_container->insert_object<PhysBasisType> (phys_space);
        }
    }

    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        // Inserting the function into the container.
        const auto obj = container->template get_const_object<Type>(id);
        full_container->insert_const_object<Type> (obj);

        // Inserting the domain of the function into the container.
        full_container->insert_const_object<DomainType> (obj->get_domain());

        // If the function is an ig function, its physical space basis is also inserted.
        const auto ig_func = dynamic_pointer_cast<const IgFunctionType>(obj);
        if (ig_func != nullptr)
            full_container->insert_const_object<PhysBasisType> (ig_func->get_basis());
    }
  });


  // Filling all physical space bases.
  PhysSpacePtrs valid_ps_ptr_types;
  for_each(valid_ps_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using DomainType = typename Type::PhysDomain;
    using RefBasisType = typename Type::RefBasis;

    // Adding const objects.
    for (const auto &id : container->template get_object_ids<Type>())
    {
        // Inserting the physical space basis into the container.
        const auto obj = container->template get_object<Type>(id);
        full_container->insert_object<Type> (obj);

        // Inserting the domain of the physical space basis into the container.
        const auto domain = const_pointer_cast<DomainType>(obj->get_physical_domain());
        Assert (domain != nullptr, ExcNullPtr());
        full_container->insert_object<DomainType> (domain);

        // Inserting the reference space basis of the physical space basis into the container.
        const auto ref_space = const_pointer_cast<RefBasisType>(obj->get_reference_basis());
        Assert (ref_space != nullptr, ExcNullPtr());
        full_container->insert_object<RefBasisType> (ref_space);
    }

    // Adding non-const objects.
    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        // Inserting the physical space basis into the container.
        const auto obj = container->template get_const_object<Type>(id);
        full_container->insert_const_object<Type> (obj);

        // Inserting the domain of the physical space basis into the container.
        full_container->insert_const_object<DomainType> (obj->get_physical_domain());

        // Inserting the reference space basis of the physical space basis into the container.
        full_container->insert_const_object<RefBasisType> (obj->get_reference_basis());
    }
  });


  // Filling all domains.
  DomainPtrs valid_dm_ptr_types;
  for_each(valid_dm_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridFuncType = typename Type::GridFuncType;

    // Adding const objects.
    for (const auto &id : container->template get_object_ids<Type>())
    {
        // Inserting the domain into the container.
        const auto obj = container->template get_object<Type>(id);
        full_container->insert_object<Type> (obj);

        // Inserting the grid function of the domain into the container.
        const auto grid_func = const_pointer_cast<GridFuncType>(obj->get_grid_function());
        Assert (grid_func != nullptr, ExcNullPtr());
        full_container->insert_object<GridFuncType> (grid_func);
    }

    // Adding non-const objects.
    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        // Inserting the domain into the container.
        const auto obj = container->template get_const_object<Type>(id);
        full_container->insert_const_object<Type> (obj);

        // Inserting the grid function of the domain into the container.
        full_container->insert_const_object<GridFuncType> (obj->get_grid_function());
    }
  });


  // Filling all grid functions.
  GridFuncPtrs valid_gr_f_ptr_types;
  for_each(valid_gr_f_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridType = typename Type::GridType;
    using IgGridFuncType = IgGridFunction<Type::dim, Type::space_dim>;
    using RefBasisType = typename IgGridFuncType::RefBasis;

    // Adding const objects.
    for (const auto &id : container->template get_object_ids<Type>())
    {
        // Inserting grid function into the container.
        const auto obj = container->template get_object<Type>(id);
        full_container->insert_object<Type> (obj);

        // Inserting the grid of the grid function into the container.
        const auto grid = const_pointer_cast<GridType>(obj->get_grid());
        full_container->insert_object<GridType> (grid);

        // If the grid function is an ig grid function, its
        // reference space basis is also inserted.
        const auto ig_g_f = dynamic_pointer_cast<IgGridFuncType>(obj);
        if (ig_g_f != nullptr)
        {
            const auto ref_space = const_pointer_cast<RefBasisType> (ig_g_f->get_basis());
            Assert (ref_space != nullptr, ExcNullPtr());
            full_container->insert_object<RefBasisType> (ref_space);
        }
    }

    // Adding non-const objects.
    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        // Inserting grid function into the container.
        const auto obj = container->template get_const_object<Type>(id);
        full_container->insert_const_object<Type> (obj);

        // Inserting the grid of the grid function into the container.
        full_container->insert_const_object<GridType> (obj->get_grid());

        // If the grid function is an ig grid function, its
        // reference space basis is also inserted.
        const auto ig_g_f = dynamic_pointer_cast<const IgGridFuncType>(obj);
        if (ig_g_f != nullptr)
            full_container->insert_const_object<RefBasisType> (ig_g_f->get_basis());
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
    using GridFuncType = GridFunction<WeightFuncType::dim, WeightFuncType::space_dim>;

    // Adding const objects.
    for (const auto &id : container->template get_object_ids<Type>())
    {
        // Inserting reference space basis function into the container.
        const auto obj = container->template get_object<Type>(id);
        full_container->insert_object<Type> (obj);

        const auto nr = dynamic_pointer_cast<NURBSType>(obj);
        const auto bs = dynamic_pointer_cast<BSplineType>(obj);

        Assert (bs != nullptr && nr != nullptr,
                ExcMessage("Invalid reference space basis type."));

        // If the reference space is a BSpline, the spline space is also inserted.
        if (bs != nullptr) // BSpline
        {
            // Inserting spline space into the container.
            const auto ss = const_pointer_cast<SpSpaceType> (bs->get_spline_space());
            full_container->insert_object<SpSpaceType> (ss);
        }
        // If the reference space is a NURBS, its BSpline space, the
        // spline space, the weight grid function, the reference space of
        // the weights and spline space space of the weights are also
        // inserted.
        else // NURBS
        {
            // Inserting spline space into the container.
            const auto ss = const_pointer_cast<SpSpaceType> (nr->get_spline_space());
            full_container->insert_object<SpSpaceType> (ss);

            // Inserting the BSpline space of the NURBS into the container.
            const auto nr_bs = const_pointer_cast<BSplineType> (nr->get_bspline_basis());
            full_container->insert_object<Type> (nr_bs);

            // Adding the weight related quantities: grid function,
            // reference space basis and its spline space.

            // Inserting the grid function of the weights.
            const auto wf = const_pointer_cast<WeightFuncType>(nr->get_weight_func());
            Assert (wf != nullptr, ExcNullPtr());
            full_container->insert_object<GridFuncType>(wf);

            // Inserting the reference space basis of the weights.
            const auto rs = const_pointer_cast<WeightFuncBasisType>(wf->get_basis());
            Assert (rs != nullptr, ExcNullPtr());
            full_container->insert_object<WeightFuncBasisType>(rs);

            // Inserting the spline space basis of the weights.
            const auto rs_ss = const_pointer_cast<WeightFuncSpSpaceType>(rs->get_spline_space());
            full_container->insert_object<WeightFuncSpSpaceType>(rs_ss);
        }
    }

    // Adding non-const objects.
    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        // Inserting reference space basis function into the container.
        const auto obj = container->template get_const_object<Type>(id);
        full_container->insert_const_object<Type> (obj);

        const auto nr = dynamic_pointer_cast<const NURBSType>(obj);
        const auto bs = dynamic_pointer_cast<const BSplineType>(obj);

        Assert (bs != nullptr && nr != nullptr,
                ExcMessage("Invalid reference space basis type."));

        // If the reference space is a BSpline, the spline space is also inserted.
        if (bs != nullptr) // BSpline
        {
            // Inserting spline space into the container.
            full_container->insert_const_object<SpSpaceType> (bs->get_spline_space());
        }
        // If the reference space is a NURBS, its BSpline space, the
        // spline space, the weight grid function, the reference space of
        // the weights and spline space space of the weights are also
        // inserted.
        else // NURBS
        {
            // Inserting spline space into the container.
            full_container->insert_const_object<SpSpaceType> (nr->get_spline_space());

            // Inserting the BSpline space of the NURBS into the container.
            full_container->insert_const_object<Type> (nr->get_bspline_basis());

            // Adding the weight related quantities: grid function,
            // reference space basis and its spline space.

            const auto wf = nr->get_weight_func();
            const auto rs = wf->get_basis();
            const auto gf = dynamic_pointer_cast<const GridFuncType> (wf);
            Assert (gf != nullptr, ExcNullPtr());

            // Inserting the grid function of the weights.
            full_container->insert_const_object<GridFuncType>(gf);

            // Inserting the reference space basis of the weights.
            full_container->insert_const_object<WeightFuncBasisType>(rs);

            // Inserting the spline space basis of the weights.
            full_container->insert_const_object<WeightFuncSpSpaceType>(rs->get_spline_space());
        }
    }
  });


  // Filling all spline spaces.
  SpSpacePtrs valid_ss_ptr_types;
  for_each(valid_ss_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridType = Grid<Type::dim>;

    // Adding const objects.
    for (const auto &id : container->template get_object_ids<Type>())
    {
        // Inserting spline space into the container.
        const auto obj = container->template get_object<Type>(id);
        full_container->insert_object<Type> (obj);

        // Inserting the grid of the spline space into the container.
        const auto grid = const_pointer_cast<GridType>(obj->get_grid());
        full_container->insert_object<GridType> (grid);
    }

    // Adding non-const objects.
    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        // Inserting spline space into the container.
        const auto obj = container->template get_const_object<Type>(id);
        full_container->insert_const_object<Type> (obj);

        // Inserting the grid of the spline space into the container.
        full_container->insert_const_object<GridType> (obj->get_grid());
    }
  });

  // Filling all grids.
  GridPtrs valid_gr_ptr_types;
  for_each(valid_gr_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    // Adding const objects.
    for (const auto &id : container->template get_object_ids<Type>())
    {
        // Inserting grid into the container.
        const auto obj = container->template get_object<Type>(id);
        full_container->insert_object<Type> (obj);
    }

    // Adding non-const objects.
    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        // Inserting grid into the container.
        const auto obj = container->template get_const_object<Type>(id);
        full_container->insert_const_object<Type> (obj);
    }
  });

  return full_container;
}

IGA_NAMESPACE_CLOSE

#endif // XML_IO
