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

#include <igatools/io/xml_document.h>
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
using std::remove_reference;
using std::dynamic_pointer_cast;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

void
ObjectsContainerWriter::
write(const string &file_path,
      const shared_ptr<ObjectsContainer> container)
{
    // Copying the objects container and filling it with all its dependencies.
    const auto full_container = shared_ptr<ObjectsContainer>(new ObjectsContainer(*container));
    full_container->fill_not_inserted_dependencies();

    const auto xml_doc = XMLDocument::create_void_document("Igatools");
    const auto igatools_elem = xml_doc->get_document_element();
    igatools_elem->add_attribute(string("FormatVersion"), string("2.0"));
    const auto new_elem = xml_doc->create_new_element("joder");
    igatools_elem->append_child_element(new_elem);

    Self_::write_grids(full_container, xml_doc);
    Self_::write_spline_spaces(full_container, xml_doc);
    Self_::write_reference_space_bases(full_container, xml_doc);
    Self_::write_grid_functions(full_container, xml_doc);
    Self_::write_domains(full_container, xml_doc);
    Self_::write_physical_space_bases(full_container, xml_doc);
    Self_::write_functions(full_container, xml_doc);

    xml_doc->write_to_file (file_path);
}



void
ObjectsContainerWriter::
write_grids (const shared_ptr<ObjectsContainer> container,
             const XMLDocPtr_ xml_doc)
{
  using GridPtrs = typename ObjectsContainer::GridPtrs;

  GridPtrs valid_grid_ptr_types;
  for_each(valid_grid_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    for (const auto &id : container->template get_object_ids<Type>())
        Self_::write_grid<Type>(container->template get_object<Type>(id),
                                xml_doc);

    for (const auto &id : container->template get_const_object_ids<Type>())
        Self_::write_grid<const Type>(container->template get_const_object<Type>(id),
                                      xml_doc);
  });
}



void
ObjectsContainerWriter::
write_spline_spaces (const shared_ptr<ObjectsContainer> container,
             const XMLDocPtr_ xml_doc)
{
  using SpSpacePtrs = typename ObjectsContainer::SpSpacePtrs;

  SpSpacePtrs valid_spline_space_ptr_types;
  for_each(valid_spline_space_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    for (const auto &id : container->template get_object_ids<Type>())
        Self_::write_spline_space<Type>(container->template get_object<Type>(id),
                                        xml_doc);

    for (const auto &id : container->template get_const_object_ids<Type>())
        Self_::write_spline_space<const Type>(container->template get_const_object<Type>(id),
                                              xml_doc);
  });
}



void
ObjectsContainerWriter::
write_reference_space_bases (const shared_ptr<ObjectsContainer> container,
             const XMLDocPtr_ xml_doc)
{
  using RefSpacePtrs = typename ObjectsContainer::RefSpacePtrs;

  RefSpacePtrs valid_ref_space_ptr_types;
  for_each(valid_ref_space_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    static const int dim = Type::dim;
    static const int range = Type::range;
    static const int rank = Type::rank;

    using NURBSType = NURBS<dim, range, rank>;
    using BSplineType = typename NURBSType::BSpBasis;

    for (const auto &id : container->template get_object_ids<Type>())
    {
        const auto ref_space = container->template get_object<Type>(id);

        const auto bs_space = dynamic_pointer_cast<BSplineType>(ref_space);
        const auto nr_space = dynamic_pointer_cast<NURBSType>(ref_space);

#ifndef NDEBUG
        AssertThrow (bs_space == nullptr && nr_space == nullptr,
                ExcMessage("Invalid reference space type."));
#endif

        if (bs_space != nullptr)
            Self_::write_bspline<BSplineType>(bs_space, xml_doc);
        else
            Self_::write_bspline<NURBSType>(nr_space, xml_doc);
    }

    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        const auto ref_space = container->template get_const_object<Type>(id);

        const auto bs_space = dynamic_pointer_cast<const BSplineType>(ref_space);
        const auto nr_space = dynamic_pointer_cast<const NURBSType>(ref_space);

#ifndef NDEBUG
        AssertThrow (bs_space == nullptr && nr_space == nullptr,
                ExcMessage("Invalid reference space type."));
#endif

        if (bs_space != nullptr)
            Self_::write_bspline<const BSplineType>(bs_space, xml_doc);
        else
            Self_::write_bspline<const NURBSType>(nr_space, xml_doc);
    }
  });
}



void
ObjectsContainerWriter::
write_grid_functions (const shared_ptr<ObjectsContainer> container,
                      const XMLDocPtr_ xml_doc)
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

    for (const auto &id : container->template get_object_ids<Type>())
    {
        const auto grid_func = container->template get_object<Type>(id);

        const auto id_f = dynamic_pointer_cast<IdGridFunc>(grid_func);
        const auto li_f = dynamic_pointer_cast<LinearGridFunc>(grid_func);
        const auto ct_f = dynamic_pointer_cast<ConstantGridFunc>(grid_func);
        const auto ig_f = dynamic_pointer_cast<IgGridFunc>(grid_func);

#ifndef NDEBUG
        AssertThrow (id_f == nullptr && li_f == nullptr &&
                ct_f == nullptr && ig_f == nullptr,
                ExcMessage("Invalid grid function type."));
#endif

        if (id_f != nullptr)
            Self_::write_identity_grid_function<IdGridFunc>(id_f, xml_doc);
        else if (li_f != nullptr)
            Self_::write_linear_grid_function<LinearGridFunc>(li_f, xml_doc);
        else if (ct_f != nullptr)
            Self_::write_constant_grid_function<ConstantGridFunc>(ct_f, xml_doc);
        else
            Self_::write_ig_grid_function<IgGridFunc>(ig_f, xml_doc);
    }

    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        const auto grid_func = container->template get_const_object<Type>(id);

        const auto id_f = dynamic_pointer_cast<const IdGridFunc>(grid_func);
        const auto li_f = dynamic_pointer_cast<const LinearGridFunc>(grid_func);
        const auto ct_f = dynamic_pointer_cast<const ConstantGridFunc>(grid_func);
        const auto ig_f = dynamic_pointer_cast<const IgGridFunc>(grid_func);

#ifndef NDEBUG
        AssertThrow (id_f == nullptr && li_f == nullptr &&
                ct_f == nullptr && ig_f == nullptr,
                ExcMessage("Invalid grid function type."));
#endif

        if (id_f != nullptr)
            Self_::write_identity_grid_function<const IdGridFunc>(id_f, xml_doc);
        else if (li_f != nullptr)
            Self_::write_linear_grid_function<const LinearGridFunc>(li_f, xml_doc);
        else if (ct_f != nullptr)
            Self_::write_constant_grid_function<const ConstantGridFunc>(ct_f, xml_doc);
        else
            Self_::write_ig_grid_function<const IgGridFunc>(ig_f, xml_doc);
    }
  });
}



void
ObjectsContainerWriter::
write_domains (const shared_ptr<ObjectsContainer> container,
               const XMLDocPtr_ xml_doc)
{
  using DomainPtrs = typename ObjectsContainer::DomainPtrs;

  DomainPtrs valid_dm_ptr_types;
  for_each(valid_dm_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    for (const auto &id : container->template get_object_ids<Type>())
        Self_::write_domain<Type>(container->template get_object<Type>(id),
                                  xml_doc);

    for (const auto &id : container->template get_const_object_ids<Type>())
        Self_::write_domain<const Type>(container->template get_const_object<Type>(id),
                                        xml_doc);
  });
}



void
ObjectsContainerWriter::
write_physical_space_bases (const shared_ptr<ObjectsContainer> container,
                            const XMLDocPtr_ xml_doc)
{
  using PhysSpacePtrs = typename ObjectsContainer::PhysSpacePtrs;

  PhysSpacePtrs valid_ps_ptr_types;
  for_each(valid_ps_ptr_types, [&](const auto &ptr_type)
  {
    using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;

    for (const auto &id : container->template get_object_ids<Type>())
        Self_::write_phys_space_basis<Type>(container->template get_object<Type>(id),
                                            xml_doc);

    for (const auto &id : container->template get_const_object_ids<Type>())
        Self_::write_phys_space_basis<const Type>(container->template get_const_object<Type>(id),
                                                  xml_doc);
  });
}



void
ObjectsContainerWriter::
write_functions (const shared_ptr<ObjectsContainer> container,
                 const XMLDocPtr_ xml_doc)
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

    for (const auto &id : container->template get_object_ids<Type>())
    {
        const auto func = container->template get_object<Type>(id);

        const auto li_f = dynamic_pointer_cast<LinearFunc>(func);
        const auto ct_f = dynamic_pointer_cast<ConstantFunc>(func);
        const auto ig_f = dynamic_pointer_cast<IgFunc>(func);

#ifndef NDEBUG
        AssertThrow (li_f == nullptr && ct_f == nullptr && ig_f == nullptr,
                ExcMessage("Invalid function type."));
#endif

        if (li_f != nullptr)
            Self_::write_linear_function<LinearFunc>(li_f, xml_doc);
        else if (ct_f != nullptr)
            Self_::write_constant_function<ConstantFunc>(ct_f, xml_doc);
        else
            Self_::write_ig_function<IgFunc>(ig_f, xml_doc);
    }

    for (const auto &id : container->template get_const_object_ids<Type>())
    {
        const auto func = container->template get_const_object<Type>(id);

        const auto li_f = dynamic_pointer_cast<const LinearFunc>(func);
        const auto ct_f = dynamic_pointer_cast<const ConstantFunc>(func);
        const auto ig_f = dynamic_pointer_cast<const IgFunc>(func);

#ifndef NDEBUG
        AssertThrow (li_f == nullptr && ct_f == nullptr && ig_f == nullptr,
                ExcMessage("Invalid function type."));
#endif

        if (li_f != nullptr)
            Self_::write_linear_function<const LinearFunc>(li_f, xml_doc);
        else if (ct_f != nullptr)
            Self_::write_constant_function<const ConstantFunc>(ct_f, xml_doc);
        else
            Self_::write_ig_function<const IgFunc>(ig_f, xml_doc);
    }
  });
}



template <class Grid>
void
ObjectsContainerWriter::
write_grid (const shared_ptr<Grid> grid,
            const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class SpSpace>
void
ObjectsContainerWriter::
write_spline_space (const shared_ptr<SpSpace> spline_space,
                    const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class BSpline>
void
ObjectsContainerWriter::
write_bspline (const shared_ptr<BSpline> bspline,
               const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class NURBS>
void
ObjectsContainerWriter::
write_nurbs (const shared_ptr<NURBS> nurbs,
             const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class IdGridFunc>
void
ObjectsContainerWriter::
write_identity_grid_function (const shared_ptr<IdGridFunc> id_func,
                              const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class ConstGridFunc>
void
ObjectsContainerWriter::
write_constant_grid_function (const shared_ptr<ConstGridFunc> const_func,
                              const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class LinearGridFunc>
void
ObjectsContainerWriter::
write_linear_grid_function (const shared_ptr<LinearGridFunc> linear_func,
                            const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class IgGridFunc>
void
ObjectsContainerWriter::
write_ig_grid_function (const shared_ptr<IgGridFunc> ig_func,
                        const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class Domain>
void
ObjectsContainerWriter::
write_domain (const shared_ptr<Domain> domain,
              const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class PhysSpaceBasis>
void
ObjectsContainerWriter::
write_phys_space_basis (const shared_ptr<PhysSpaceBasis> phys_space,
                        const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class IgFunction>
void
ObjectsContainerWriter::
write_ig_function (const shared_ptr<IgFunction> ig_function,
                   const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class ConstantFunction>
void
ObjectsContainerWriter::
write_constant_function (const shared_ptr<ConstantFunction> const_function,
                         const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}



template <class LinearFunction>
void
ObjectsContainerWriter::
write_linear_function (const shared_ptr<LinearFunction> linear_function,
                       const XMLDocPtr_ xml_doc)
{
    AssertThrow (false, ExcNotImplemented());
}

IGA_NAMESPACE_CLOSE

#endif // XML_IO
