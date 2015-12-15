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

#include <igatools/io/objects_container_xml_reader.h>

#ifdef XML_IO

#include <igatools/io/objects_container_reader-XML_schema.h>

#include <igatools/io/xml_document.h>
#include <igatools/io/xml_element.h>
#include <igatools/base/objects_container.h>

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
using std::to_string;
using std::shared_ptr;
using std::set;
using std::remove_reference;
using std::dynamic_pointer_cast;
using std::copy;
using std::inserter;

IGA_NAMESPACE_OPEN

shared_ptr<ObjectsContainer>
ObjectsContainerXMLReader::
parse(const string &file_path)
{
  const auto xml_doc = XMLDocument::parse_from_file(file_path,
                       ObjectsContainerXMLReader::XML_SCHEMA_);
  const auto xml_elem = xml_doc->get_document_element();

  const auto container = ObjectsContainer::create();

  // Checking for repeated iga object ids.
  SafeSTLSet<Index> object_ids;
  for (const auto &el : xml_elem->get_children_elements())
  {
    const Index obj_id = el->get_attribute<Index>("LocalObjectId");
    AssertThrow(object_ids.find(obj_id) == object_ids.cend(),
                ExcMessage("LocalObjectId " + to_string(obj_id) + " is "
                           "used more than once."));
    object_ids.insert(obj_id);
  }
  IdMap_ id_map;

  const bool parse_as_constant = false;
  Self_::parse_grids(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_spline_spaces(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_bsplines(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_grid_functions_and_nurbs(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_domains(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_phys_spaces(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_functions(xml_elem, parse_as_constant, id_map, container);

  return container;
}



shared_ptr<ObjectsContainer>
ObjectsContainerXMLReader::
parse_const(const string &file_path)
{
  const auto xml_doc = XMLDocument::parse_from_file(file_path,
                       ObjectsContainerXMLReader::XML_SCHEMA_);
  const auto xml_elem = xml_doc->get_document_element();

  const auto container = ObjectsContainer::create();

  // Checking for repeated iga object ids.
  SafeSTLSet<Index> object_ids;
  for (const auto &el : xml_elem->get_children_elements())
  {
    const Index obj_id = el->get_attribute<Index>("LocalObjectId");
    AssertThrow(object_ids.find(obj_id) == object_ids.cend(),
                ExcMessage("LocalObjectId " + to_string(obj_id) + " is "
                           "used more than once."));
    object_ids.insert(obj_id);
  }
  IdMap_ id_map;

  const bool parse_as_constant = true;
  Self_::parse_grids(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_spline_spaces(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_bsplines(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_grid_functions_and_nurbs(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_domains(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_phys_spaces(xml_elem, parse_as_constant, id_map, container);
  Self_::parse_functions(xml_elem, parse_as_constant, id_map, container);

  return container;
}



void
ObjectsContainerXMLReader::
parse_grids(const shared_ptr<XMLElement> xml_elem,
            const bool parse_as_constant,
            IdMap_ &id_map,
            const shared_ptr<ObjectsContainer> container)
{
  const auto grid_elems = xml_elem->get_children_elements("Grid");
  for (const auto &ge : grid_elems)
  {
    const int grid_dim = ge->get_attribute<int>("Dim");

    using GridPtrs = typename ObjectsContainer::GridPtrs;
    GridPtrs valid_grid_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_grid_ptr_types, [&](const auto &grid_ptr_type)
    {
      if (found)
        return;

      using GridType = typename remove_reference<decltype(grid_ptr_type)>::type::element_type;
      static const int dim = GridType::dim;

      if (grid_dim == dim)
      {
        found = true;
        parse_grid<dim>(ge, parse_as_constant, id_map, container);
      }
    });

    // Grid dimension not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("Grid",
                                                     ge->get_attribute<Index>("LocalObjectId"),
                                                     SafeSTLVector<int>(1, grid_dim))
                           + " is not a valid type. Possibly the type was not "
                           "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_spline_spaces(const shared_ptr<XMLElement> xml_elem,
                    const bool parse_as_constant,
                    IdMap_ &id_map,
                    const shared_ptr<ObjectsContainer> container)
{
  for (const auto &ssp : xml_elem->get_children_elements("SplineSpace"))
  {
    const int ssp_dim = ssp->get_attribute<int>("Dim");
    const int ssp_range = ssp->get_attribute<int>("Range");
    const int ssp_rank = ssp->get_attribute<int>("Rank");

    using SpSpacePtrs = typename ObjectsContainer::SpSpacePtrs;
    SpSpacePtrs valid_rsp_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_rsp_ptr_types, [&](const auto &ssp_ptr_type)
    {
      if (found)
        return;

      using SplineSpaceType = typename
                              remove_reference<decltype(ssp_ptr_type)>::type::element_type;
      static const int dim = SplineSpaceType::dim;
      static const int range = SplineSpaceType::range;
      static const int rank = SplineSpaceType::rank;

      if (ssp_dim == dim && ssp_range == range && ssp_rank == rank)
      {
        found = true;
        parse_spline_space<dim, range, rank>(ssp, parse_as_constant, id_map, container);
      }
    });

    // SplineSpace dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("SplineSpace",
                                                     ssp->get_attribute<Index>("LocalObjectId"),
    {{ssp_dim, ssp_range, ssp_rank}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_bsplines(const shared_ptr<XMLElement> xml_elem,
               const bool parse_as_constant,
               IdMap_ &id_map,
               const shared_ptr<ObjectsContainer> container)
{
  for (const auto &bs : xml_elem->get_children_elements("BSpline"))
  {
    const int bs_dim = bs->get_attribute<int>("Dim");
    const int bs_range = bs->get_attribute<int>("Range");
    const int bs_rank = bs->get_attribute<int>("Rank");

    using BSplinePtrs = typename ObjectsContainer::RefSpacePtrs;
    BSplinePtrs valid_bs_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_bs_ptr_types, [&](const auto &bs_ptr_type)
    {
      if (found)
        return;

      using BSplineType = typename
                          remove_reference<decltype(bs_ptr_type)>::type::element_type;
      static const int dim = BSplineType::dim;
      static const int range = BSplineType::range;
      static const int rank = BSplineType::rank;

      if (bs_dim == dim && bs_range == range && bs_rank == rank)
      {
        found = true;
        parse_bspline<dim, range, rank>(bs, parse_as_constant, id_map, container);
      }
    });

    // BSpline dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("BSpline",
                                                     bs->get_attribute<Index>("LocalObjectId"),
    {{bs_dim, bs_range, bs_rank}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_nurbs(const shared_ptr<XMLElement> xml_elem,
            const bool parse_as_constant,
            IdMap_ &id_map,
            const shared_ptr<ObjectsContainer> container)
{
  for (const auto &nr : xml_elem->get_children_elements("NURBS"))
  {
    const int nr_dim = nr->get_attribute<int>("Dim");
    const int nr_range = nr->get_attribute<int>("Range");
    const int nr_rank = nr->get_attribute<int>("Rank");

    using NURBSPtrs = typename ObjectsContainer::RefSpacePtrs;
    NURBSPtrs valid_nr_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_nr_ptr_types, [&](const auto &nr_ptr_type)
    {
      if (found)
        return;

      using NURBSType = typename
                        remove_reference<decltype(nr_ptr_type)>::type::element_type;
      static const int dim   = NURBSType::dim;
      static const int range = NURBSType::range;
      static const int rank  = NURBSType::rank;

      if (nr_dim == dim && nr_range == range && nr_rank == rank)
      {
        found = true;
        parse_nurbs<dim, range, rank>(nr, parse_as_constant, id_map, container);
      }
    });

    // NURBS dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("NURBS",
                                                     nr->get_attribute<Index>("LocalObjectId"),
    {{nr_dim, nr_range, nr_rank}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_identity_grid_functions(const shared_ptr<XMLElement> xml_elem,
                              const bool parse_as_constant,
                              IdMap_ &id_map,
                              const shared_ptr<ObjectsContainer> container)
{
  for (const auto &id : xml_elem->get_children_elements("IdentityGridFunction"))
  {
    const int id_dim = id->get_attribute<int>("Dim");

    using GridPtrs = typename ObjectsContainer::GridPtrs;
    GridPtrs valid_id_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_id_ptr_types, [&](const auto &id_ptr_type)
    {
      if (found)
        return;

      using GridType = typename
                       remove_reference<decltype(id_ptr_type)>::type::element_type;
      static const int dim   = GridType::dim;

      if (id_dim == dim)
      {
        found = true;
        parse_identity_grid_function<dim>(id, parse_as_constant, id_map, container);
      }
    });

    // IdentityGridFunction dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("IdentityGridFunction",
                                                     id->get_attribute<Index>("LocalObjectId"),
                                                     SafeSTLVector<int>(1, id_dim))
                           + " is not a valid type. Possibly the type was not "
                           "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_constant_grid_functions(const shared_ptr<XMLElement> xml_elem,
                              const bool parse_as_constant,
                              IdMap_ &id_map,
                              const shared_ptr<ObjectsContainer> container)
{
  for (const auto &cgf : xml_elem->get_children_elements("ConstantGridFunction"))
  {
    const int cgf_dim = cgf->get_attribute<int>("Dim");
    const int cgf_space_dim = cgf->get_attribute<int>("Spacedim"); // This is going to change.

    using GridFunctionPtrs = typename ObjectsContainer::GridFuncPtrs;
    GridFunctionPtrs valid_cgf_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_cgf_ptr_types, [&](const auto &cgf_ptr_type)
    {
      if (found)
        return;

      using GridFuncType = typename
                           remove_reference<decltype(cgf_ptr_type)>::type::element_type;
      static const int dim   = GridFuncType::dim;
      static const int space_dim  = GridFuncType::space_dim;

      if (cgf_dim == dim && cgf_space_dim == space_dim)
      {
        found = true;
        parse_constant_grid_function<dim, space_dim>(cgf, parse_as_constant, id_map, container);
      }
    });

    // ConstantGridFunction dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("ConstantGridFunction",
                                                     cgf->get_attribute<Index>("LocalObjectId"),
    {{cgf_dim, cgf_space_dim}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_linear_grid_functions(const shared_ptr<XMLElement> xml_elem,
                            const bool parse_as_constant,
                            IdMap_ &id_map,
                            const shared_ptr<ObjectsContainer> container)
{
  for (const auto &cgf : xml_elem->get_children_elements("LinearGridFunction"))
  {
    const int cgf_dim = cgf->get_attribute<int>("Dim");
    const int cgf_space_dim = cgf->get_attribute<int>("Spacedim"); // This is going to change.

    using GridFunctionPtrs = typename ObjectsContainer::GridFuncPtrs;
    GridFunctionPtrs valid_cgf_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_cgf_ptr_types, [&](const auto &cgf_ptr_type)
    {
      if (found)
        return;

      using GridFuncType = typename
                           remove_reference<decltype(cgf_ptr_type)>::type::element_type;
      static const int dim   = GridFuncType::dim;
      static const int space_dim  = GridFuncType::space_dim;

      if (cgf_dim == dim && cgf_space_dim == space_dim)
      {
        found = true;
        parse_linear_grid_function<dim, space_dim>(cgf, parse_as_constant, id_map, container);
      }
    });

    // LinearGridFunction dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("LinearGridFunction",
                                                     cgf->get_attribute<Index>("LocalObjectId"),
    {{cgf_dim, cgf_space_dim}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_ig_grid_functions(const shared_ptr<XMLElement> xml_elem,
                        const bool parse_as_constant,
                        const bool &first_parsing,
                        IdMap_ &id_map,
                        const shared_ptr<ObjectsContainer> container)
{
  for (const auto &gf : xml_elem->get_children_elements("IgGridFunction"))
  {
    const int gf_dim = gf->get_attribute<int>("Dim");
    const int gf_space_dim = gf->get_attribute<int>("Spacedim"); // This is going to change.

    using GridFunctionPtrs = typename ObjectsContainer::GridFuncPtrs;
    GridFunctionPtrs valid_gf_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_gf_ptr_types, [&](const auto &gf_ptr_type)
    {
      if (found)
        return;

      using GridFuncType = typename
                           remove_reference<decltype(gf_ptr_type)>::type::element_type;
      static const int dim   = GridFuncType::dim;
      static const int space_dim   = GridFuncType::space_dim;

      if (gf_dim == dim && gf_space_dim == space_dim)
      {
        found = true;
        parse_ig_grid_function<dim, space_dim>(gf, parse_as_constant, first_parsing, id_map, container);
      }
    });

    // NURBS dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("IgGridFunction",
                                                     gf->get_attribute<Index>("LocalObjectId"),
    {{gf_dim, gf_space_dim}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_grid_functions_and_nurbs(const shared_ptr<XMLElement> xml_elem,
                               const bool parse_as_constant,
                               IdMap_ &id_map,
                               const shared_ptr<ObjectsContainer> container)
{
  // Due to the relationship between ig grid functions and NURBS,
  // their parsing must be done in a specific order. That is:
  //
  //   1 - parse grid functions that are not ig.
  //   2 - parse ig grid functions built upon a BSpline.
  //   3 - parse nurbs
  //   4 - parse the remaning ig grid functions.


  // Parsing identity grid functions.
  parse_identity_grid_functions(xml_elem, parse_as_constant, id_map, container);

  // Parsing constant grid functions.
  parse_constant_grid_functions(xml_elem, parse_as_constant, id_map, container);

  // Parsing linear grid functions.
  parse_linear_grid_functions(xml_elem, parse_as_constant, id_map, container);

  // Parsing ig grid functions built upon a BSpline.
  bool first_parsing = true;
  parse_ig_grid_functions(xml_elem, parse_as_constant, first_parsing, id_map, container);

  // Parsing NURBS.
  Self_::parse_nurbs(xml_elem, parse_as_constant, id_map, container);

  // Parsing the remaining ig grid functions.
  first_parsing = false;
  parse_ig_grid_functions(xml_elem, parse_as_constant, first_parsing, id_map, container);
}



void
ObjectsContainerXMLReader::
parse_domains(const shared_ptr<XMLElement> xml_elem,
              const bool parse_as_constant,
              IdMap_ &id_map,
              const shared_ptr<ObjectsContainer> container)
{
  for (const auto &dm : xml_elem->get_children_elements("Domain"))
  {
    const int dm_dim = dm->get_attribute<int>("Dim");
    const int dm_codim = dm->get_attribute<int>("Codim");

    using DomainPtrs = typename ObjectsContainer::DomainPtrs;
    DomainPtrs valid_dm_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_dm_ptr_types, [&](const auto &dm_ptr_type)
    {
      if (found)
        return;

      using DomainType = typename
                         remove_reference<decltype(dm_ptr_type)>::type::element_type;
      static const int dim = DomainType::dim;
      static const int space_dim = DomainType::space_dim;
      static const int codim = space_dim - dim;

      if (dm_dim == dim && dm_codim == codim)
      {
        found = true;
        parse_domain<dim, codim>(dm, parse_as_constant, id_map, container);
      }
    });

    // Domains dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("Domain",
                                                     dm->get_attribute<Index>("LocalObjectId"),
    {{dm_dim, dm_codim}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_phys_spaces(const shared_ptr<XMLElement> xml_elem,
                  const bool parse_as_constant,
                  IdMap_ &id_map,
                  const shared_ptr<ObjectsContainer> container)
{
  for (const auto &ps : xml_elem->get_children_elements("PhysicalSpaceBasis"))
  {
    const int ps_dim = ps->get_attribute<int>("Dim");
    const int ps_codim = ps->get_attribute<int>("Codim");
    const int ps_range = ps->get_attribute<int>("Range");
    const int ps_rank = ps->get_attribute<int>("Rank");

    using PhysSpacePtrs = typename ObjectsContainer::PhysSpacePtrs;
    PhysSpacePtrs valid_ps_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_ps_ptr_types, [&](const auto &ps_ptr_type)
    {
      if (found)
        return;

      using PhysSpaceType = typename
                            remove_reference<decltype(ps_ptr_type)>::type::element_type;
      static const int dim = PhysSpaceType::dim;
      static const int range = PhysSpaceType::range;
      static const int rank = PhysSpaceType::rank;
      static const int codim = PhysSpaceType::codim;

      if (ps_dim == dim && ps_range == range && ps_rank == rank && ps_codim == codim)
      {
        found = true;
        parse_phys_space<dim, codim, range, rank>(ps, parse_as_constant, id_map, container);
      }
    });

    // PhysicalSpaceBasis dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("PhysicalSpaceBasis",
                                                     ps->get_attribute<Index>("LocalObjectId"),
    {{ps_dim, ps_range, ps_rank, ps_codim}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_functions(const shared_ptr<XMLElement> xml_elem,
                const bool parse_as_constant,
                IdMap_ &id_map,
                const shared_ptr<ObjectsContainer> container)
{
  parse_ig_functions(xml_elem, parse_as_constant, id_map, container);
  parse_constant_functions(xml_elem, parse_as_constant, id_map, container);
  parse_linear_functions(xml_elem, parse_as_constant, id_map, container);
}



void
ObjectsContainerXMLReader::
parse_ig_functions(const shared_ptr<XMLElement> xml_elem,
                   const bool parse_as_constant,
                   IdMap_ &id_map,
                   const shared_ptr<ObjectsContainer> container)
{
  for (const auto &fn : xml_elem->get_children_elements("IgFunction"))
  {
    const int fn_dim = fn->get_attribute<int>("Dim");
    const int fn_codim = fn->get_attribute<int>("Codim");
    const int fn_range = fn->get_attribute<int>("Range");
    const int fn_rank = fn->get_attribute<int>("Rank");

    using FunctionPtrs = typename ObjectsContainer::FunctionPtrs;
    FunctionPtrs valid_fn_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_fn_ptr_types, [&](const auto &fn_ptr_type)
    {
      if (found)
        return;

      using FunctionType = typename
                           remove_reference<decltype(fn_ptr_type)>::type::element_type;
      static const int dim = FunctionType::dim;
      static const int range = FunctionType::range;
      static const int rank = FunctionType::rank;
      static const int codim = FunctionType::codim;

      if (fn_dim == dim && fn_range == range && fn_rank == rank && fn_codim == codim)
      {
        found = true;
        parse_ig_function<dim, codim, range, rank>(fn, parse_as_constant, id_map, container);
      }
    });

    // Function dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("IgFunction",
                                                     fn->get_attribute<Index>("LocalObjectId"),
    {{fn_dim, fn_codim, fn_range, fn_rank}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_constant_functions(const shared_ptr<XMLElement> xml_elem,
                         const bool parse_as_constant,
                         IdMap_ &id_map,
                         const shared_ptr<ObjectsContainer> container)
{
  for (const auto &fn : xml_elem->get_children_elements("ConstantFunction"))
  {
    const int fn_dim = fn->get_attribute<int>("Dim");
    const int fn_codim = fn->get_attribute<int>("Codim");
    const int fn_range = fn->get_attribute<int>("Range");
    const int fn_rank = fn->get_attribute<int>("Rank");

    using FunctionPtrs = typename ObjectsContainer::FunctionPtrs;
    FunctionPtrs valid_fn_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_fn_ptr_types, [&](const auto &fn_ptr_type)
    {
      if (found)
        return;

      using FunctionType = typename
                           remove_reference<decltype(fn_ptr_type)>::type::element_type;
      static const int dim = FunctionType::dim;
      static const int range = FunctionType::range;
      static const int rank = FunctionType::rank;
      static const int codim = FunctionType::codim;

      if (fn_dim == dim && fn_range == range && fn_rank == rank && fn_codim == codim)
      {
        found = true;
        parse_constant_function<dim, codim, range, rank>(fn, parse_as_constant, id_map, container);
      }
    });

    // Function dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("ConstantFunction",
                                                     fn->get_attribute<Index>("LocalObjectId"),
    {{fn_dim, fn_codim, fn_range, fn_rank}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



void
ObjectsContainerXMLReader::
parse_linear_functions(const shared_ptr<XMLElement> xml_elem,
                       const bool parse_as_constant,
                       IdMap_ &id_map,
                       const shared_ptr<ObjectsContainer> container)
{
  for (const auto &fn : xml_elem->get_children_elements("LinearFunction"))
  {
    const int fn_dim = fn->get_attribute<int>("Dim");
    const int fn_codim = fn->get_attribute<int>("Codim");
    const int fn_range = fn->get_attribute<int>("Range");

    using FunctionPtrs = typename ObjectsContainer::FunctionPtrs;
    FunctionPtrs valid_fn_ptr_types;

    bool found = false;
    boost::fusion::for_each(valid_fn_ptr_types, [&](const auto &fn_ptr_type)
    {
      if (found)
        return;

      using FunctionType = typename
                           remove_reference<decltype(fn_ptr_type)>::type::element_type;
      static const int dim = FunctionType::dim;
      static const int range = FunctionType::range;
      static const int codim = FunctionType::codim;

      if (fn_dim == dim && fn_range == range && fn_codim == codim)
      {
        found = true;
        parse_linear_function<dim, codim, range>(fn, parse_as_constant, id_map, container);
      }
    });

    // Function dimensions not found
    AssertThrow(found,
                ExcMessage(Self_::get_type_id_string("LinearFunction",
                                                     fn->get_attribute<Index>("LocalObjectId"),
    {{fn_dim, fn_codim, fn_range}})
    + " is not a valid type. Possibly the type was not "
    "instantiated for the specified dimensions."));
  }
}



template <int dim>
void
ObjectsContainerXMLReader::
parse_grid(const shared_ptr<XMLElement> xml_elem,
           const bool parse_as_constant,
           IdMap_ &id_map,
           const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "Grid", ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));

  using GridType = Grid<dim>;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("Grid", local_object_id,
                                                       SafeSTLVector<int>(1, dim));

  const auto knots_children = xml_elem->get_children_elements("Knots");
  // Checking the number of knot vector match with the dimension.
  AssertThrow(dim == knots_children.size(),
              ExcMessage("Parsing " + parsing_msg + ", the number of "
                         "knots vectors is not valid."));

  SafeSTLArray<SafeSTLVector<Real>, dim> knots;

  SafeSTLSet<Index> parsed_dirs;

  for (const auto &ke : knots_children)
  {
    const auto dir = ke->get_attribute<Index>("Direction");
    // Checking the <= 0 direction < dim
    AssertThrow(dir >= 0 && dir < dim,
                ExcMessage("Parsing knot vectors for " + parsing_msg +
                           ", not valid Direction=" + to_string(dir) + "."));
    // Checking the direction has not been defined before.
    AssertThrow(parsed_dirs.find(dir) == parsed_dirs.cend(),
                ExcMessage("Parsing knot vectors for " + parsing_msg +
                           ", Direction=" + to_string(dir) + " defined"
                           " more than once."));
    parsed_dirs.insert(dir);

    knots[dir] = ke->get_values_vector<Real>();

    const auto size = ke->get_attribute<int>("Size");
    // Checking that the specified size matches with the actual vector size.
    AssertThrow(size == knots[dir].size(),
                ExcMessage("Parsing knot vectors for " + parsing_msg +
                           ", in Direction=" + to_string(dir) +
                           " Size=" + to_string(size) + " do not match "
                           "with the vector size."));
  }

  if (parse_as_constant)
  {
    const auto grid = GridType::const_create(knots);
    const auto unique_id = grid->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<GridType>(grid);
  }
  else
  {
    const auto grid = GridType::create(knots);
    const auto unique_id = grid->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<GridType>(grid);

  }
}



template <int dim, int range, int rank>
void
ObjectsContainerXMLReader::
parse_spline_space(const shared_ptr<XMLElement> xml_elem,
                   const bool parse_as_constant,
                   IdMap_ &id_map,
                   const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "SplineSpace",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Range") == range,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
  Assert(xml_elem->get_attribute<int>("Rank") == rank,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));

  using GridType = Grid<dim>;
  using SplineSpaceType = SplineSpace<dim, range, rank>;
  using DegreeTable = typename SplineSpaceType::DegreeTable;
  using MultiplicityTable = typename SplineSpaceType::MultiplicityTable;
  using PeriodicityTable = typename SplineSpaceType::PeriodicityTable;
  static const int n_components = SplineSpaceType::n_components;

  DegreeTable deg_table;
  MultiplicityTable mult_table;

  PeriodicityTable period_table;
  // Initializing default periodicity
  const bool default_periodicity = false;
  for (auto &p_t : period_table)
    for (auto &p : p_t)
      p = default_periodicity;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const auto grid_tag = xml_elem->get_single_element("Grid");
  const auto local_grid_id = grid_tag->get_attribute<Index>("GetFromLocalObjectId");

  const string parsing_msg = Self_::get_type_id_string("SplineSpace", local_object_id,
  {{dim, range, rank}});

  // Checking the grid with proper dimension and id exists.
  AssertThrow(id_map.find(local_grid_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<GridType> (id_map.at(local_grid_id)) :
               container->is_object_present<GridType> (id_map.at(local_grid_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching definition" +
                         "for " + Self_::get_type_id_string("Grid", local_grid_id,
                                                            SafeSTLVector<Index>(dim)) + "."));

  SharedPtrConstnessHandler<GridType> grid;
  if (parse_as_constant)
    grid = SharedPtrConstnessHandler<GridType>
           (container->get_const_object<GridType>(id_map.at(local_grid_id)));
  else
    grid = SharedPtrConstnessHandler<GridType>
           (container->get_object<GridType>(id_map.at(local_grid_id)));

  const auto grid_num_intervals = grid->get_num_intervals();

  // Parsing spline space components.
  const auto comps_elem = xml_elem->get_single_element("SplineSpaceComponents");

  // Checking that the right number of components is defined.
  const auto comp_children = comps_elem->get_children_elements("SplineSpaceComponent");
  AssertThrow(n_components == comp_children.size(),
              ExcMessage("Parsing " + parsing_msg + ", the number of "
                         "SplineSpaceComponent XML elements does not match "
                         "with the number of components of the space."));

  SafeSTLSet<Index> parsed_comps;
  for (const auto &comp_elem : comp_children)
  {
    const auto comp_id = comp_elem->get_attribute<Index>("ComponentId");
    // Checking the <= 0 comp_id < n_components
    AssertThrow(comp_id >= 0 && comp_id < n_components,
                ExcMessage("Parsing SplineSpaceComponents for " + parsing_msg +
                           ", not valid ComponentId=" + to_string(comp_id) + "."));
    // Checking the component id has not been defined before.
    AssertThrow(parsed_comps.find(comp_id) == parsed_comps.cend(),
                ExcMessage("Parsing SplineSpaceComponents for " + parsing_msg +
                           ", ComponentId=" + to_string(comp_id) + " defined"
                           " more than once."));
    parsed_comps.insert(comp_id);

    const auto degree_elem = comp_elem->get_single_element("Degrees");
    const auto degs_vector = degree_elem->get_values_vector<Index>();
    // Check here that degs_vector.size() == dim
    AssertThrow(degs_vector.size() == dim,
                ExcMessage("Parsing " + parsing_msg +
                           ", in SplineSpaceComponent ComponentId=" +
                           to_string(comp_id) + " the number of "
                           "degrees does not match with the Space dimension."));
    for (int d = 0; d < dim; ++d)
      deg_table[comp_id][d] = degs_vector[d];

    const auto int_mults_elem =
      comp_elem->get_single_element("InteriorMultiplicities")
      ->get_children_elements("InteriorMultiplicities");
    // Check here that int_muls_elem.size() == dim
    AssertThrow(int_mults_elem.size() == dim,
                ExcMessage("Parsing " + parsing_msg +
                           ", in SplineSpaceComponent ComponentId=" +
                           to_string(comp_id) + " the number of "
                           "InteriorMultiplicities XML elements does not match "
                           "with the Space dimension."));

    SafeSTLSet<Index> parsed_dirs;
    for (const auto &im : int_mults_elem)
    {
      const auto dir = im->get_attribute<Index>("Direction");

      // Checking the <= 0 direction < dim
      AssertThrow(dir >= 0 && dir < dim,
                  ExcMessage("Parsing InteriorMultiplicities for " +
                             parsing_msg +
                             " SplineSpaceComponent ComponentId=" +
                             to_string(comp_id) +", not valid "
                             "Direction=" + to_string(dir) + "."));

      // Checking the direction has not been defined before.
      AssertThrow(parsed_dirs.find(dir) == parsed_dirs.cend(),
                  ExcMessage("Parsing InteriorMultiplicities for " +
                             parsing_msg +
                             " SplineSpaceComponent ComponentId=" +
                             to_string(comp_id) +", Direction=" +
                             to_string(dir) + " defined more than once."));
      parsed_dirs.insert(dir);

      const auto mults = im->get_values_vector<Index>();
      mult_table[comp_id].copy_data_direction(dir, mults);
      const auto size = im->get_attribute<int>("Size");

      // Checking that the specified size matches with the actual vector size.
      AssertThrow(size == mults.size(),
                  ExcMessage("Parsing InteriorMultiplicities for " +
                             parsing_msg +
                             " SplineSpaceComponent ComponentId=" +
                             to_string(comp_id) +" in Direction=" +
                             to_string(dir) + ", Size=" +
                             to_string(size) + " do not match "
                             "with the vector size."));

      // Check here that the multiplicities match with the grid.
      AssertThrow((grid_num_intervals[dir] - 1) == mults.size(),
                  ExcMessage("Parsing InteriorMultiplicities for " +
                             parsing_msg +
                             " SplineSpaceComponent ComponentId=" +
                             to_string(comp_id) +" in Direction=" +
                             to_string(dir) + ", the number of "
                             "multiplicities do not match with the"
                             " grid intervals."));
    }

    // Parsing periodicity, if defined.
    if (comp_elem->has_element("Periodicity"))
    {
      const auto periodic_vector = comp_elem->
                                   get_single_element("Periodicity")->get_values_vector<bool>();

      // Checking the dimensions of the periodicities vector.
      AssertThrow(periodic_vector.size() == dim,
                  ExcMessage("Parsing " + parsing_msg +
                             ", in SplineSpaceComponent ComponentId=" +
                             to_string(comp_id) + " the number of "
                             " values defined for Periodicity does not "
                             "match with the Space dimension."));

      for (int d = 0; d < dim; ++d)
        period_table[comp_id][d] = periodic_vector[d];
    }
  } // Spline Space components

  if (parse_as_constant)
  {
    const auto spline_space = SplineSpaceType::const_create
                              (deg_table, grid.get_ptr_const_data(), mult_table, period_table);
    const auto unique_id = spline_space->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<SplineSpaceType>(spline_space);
  }
  else
  {
    const auto spline_space = SplineSpaceType::create
                              (deg_table, grid.get_ptr_data(), mult_table, period_table);
    const auto unique_id = spline_space->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<SplineSpaceType>(spline_space);
  }

}



template <int dim, int range, int rank>
void
ObjectsContainerXMLReader::
parse_bspline(const shared_ptr<XMLElement> xml_elem,
              const bool parse_as_constant,
              IdMap_ &id_map,
              const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "BSpline",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Range") == range,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
  Assert(xml_elem->get_attribute<int>("Rank") == rank,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));

  using BSplineType = BSpline<dim, range, rank>;
  using RefSpaceType = ReferenceSpaceBasis<dim, range, rank>;
  using SplineSpaceType = SplineSpace<dim, range, rank>;
  using EndBehaviourTable = typename SplineSpaceType::EndBehaviourTable;
  static const int n_components = SplineSpaceType::n_components;

  EndBehaviourTable end_beh_table;
  // Initializing to default values.
  for (auto &eb_c : end_beh_table)
    for (auto &eb : eb_c)
      eb = BasisEndBehaviour::interpolatory;


  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const auto ssp_tag = xml_elem->get_single_element("SplineSpace");
  const auto local_ssp_id = ssp_tag->get_attribute<Index>("GetFromLocalObjectId");

  const string parsing_msg = Self_::get_type_id_string("BSpline", local_object_id,
  {{dim, range, rank}});

  // Checking the spline space with proper dimension and id exists.
  AssertThrow(id_map.find(local_ssp_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<SplineSpaceType> (id_map.at(local_ssp_id)) :
               container->is_object_present<SplineSpaceType> (id_map.at(local_ssp_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("SplineSpace", local_ssp_id,
  {{dim, range, rank}}) + "."));

  SharedPtrConstnessHandler<SplineSpaceType> ssp;
  if (parse_as_constant)
    ssp = SharedPtrConstnessHandler<SplineSpaceType>
          (container->get_const_object<SplineSpaceType>(id_map.at(local_ssp_id)));
  else
    ssp = SharedPtrConstnessHandler<SplineSpaceType>
          (container->get_object<SplineSpaceType>(id_map.at(local_ssp_id)));

  // Parsing end bevaviour, if exists
  if (xml_elem->has_element("EndBehaviour"))
  {

    // Checking that the right number of components is defined.
    const auto eb_elems = xml_elem->get_single_element("EndBehaviour")
                          ->get_children_elements("EndBehaviour");
    AssertThrow(n_components == eb_elems.size(),
                ExcMessage("Parsing " + parsing_msg + ", the number of "
                           "EndBehaviour elements does not match "
                           "with the number of components of the space."));

    const auto &ssp_periodic_table = ssp->get_periodic_table();

    SafeSTLSet<Index> parsed_comps;
    for (const auto &eb : eb_elems)
    {
      const auto comp_id = eb->get_attribute<Index>("ComponentId");
      // Checking the <= 0 comp_id < n_components
      AssertThrow(comp_id >= 0 && comp_id < n_components,
                  ExcMessage("Parsing EndBehaviour for " + parsing_msg +
                             ", not valid ComponentId=" + to_string(comp_id) + "."));
      // Checking the component id has not been defined before.
      AssertThrow(parsed_comps.find(comp_id) == parsed_comps.cend(),
                  ExcMessage("Parsing EndBehaviour for " + parsing_msg +
                             ", ComponentId=" + to_string(comp_id) + " defined"
                             " more than once."));
      parsed_comps.insert(comp_id);

      const auto string_vec =  eb->get_values_vector<string>();

      // Checking the dimension of the end behaviour vector.
      AssertThrow(string_vec.size() == dim,
                  ExcMessage("Parsing " + parsing_msg +
                             ", in EndBehaviour ComponentId=" +
                             to_string(comp_id) + " the number of "
                             " values defined does not match with the "
                             "Space dimension."));

      for (int dir = 0; dir < dim; ++dir)
      {
        const auto &sv = string_vec[dir];
        const bool ssp_periodic = ssp_periodic_table[comp_id][dir];
        if (sv == "interpolatory")
        {
          AssertThrow(!ssp_periodic,
                      ExcMessage("Parsing " + parsing_msg + ", in "
                                 "EndBehaviour ComponentId=" + to_string(comp_id) +
                                 " Direction=" + to_string(dir) + ", behaviour " +
                                 sv + " do not match with Spline Space, that is "
                                 "periodic for this component and direction."));
          end_beh_table[comp_id][dir] = BasisEndBehaviour::interpolatory;
        }
        else if (sv == "end_knots")
        {
          AssertThrow(false,
                      ExcMessage("Parsing " + parsing_msg + ", in "
                                 "EndBehaviour ComponentId=" + to_string(comp_id) +
                                 " Direction=" + to_string(dir) + ", behaviour " +
                                 sv + ". CURRENTLY end_knots cannot be selected from"
                                 " the XML input file."));

          AssertThrow(!ssp_periodic,
                      ExcMessage("Parsing " + parsing_msg + ", in "
                                 "EndBehaviour ComponentId=" + to_string(comp_id) +
                                 " Direction=" + to_string(dir) + ", behaviour " +
                                 sv + " do not match with Spline Space, that is "
                                 "periodic for this component and direction."));
          end_beh_table[comp_id][dir] = BasisEndBehaviour::end_knots;
        }
        else if (sv == "periodic")
        {
          AssertThrow(ssp_periodic,
                      ExcMessage("Parsing " + parsing_msg + ", in "
                                 "EndBehaviour ComponentId=" + to_string(comp_id) +
                                 " Direction=" + to_string(dir) + ", behaviour " +
                                 sv + " do not match with Spline Space, that is not "
                                 "periodic for this component and direction."));
          end_beh_table[comp_id][dir] = BasisEndBehaviour::periodic;
        }
        // If is not one of those types, the XML schema should have
        // thrown the error before.
      }
    }
  }

  if (parse_as_constant)
  {
    const auto bspline = BSplineType::const_create(ssp.get_ptr_const_data(), end_beh_table);
    const auto unique_id = bspline->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<RefSpaceType>(bspline);
  }
  else
  {
    const auto bspline = BSplineType::create(ssp.get_ptr_data(), end_beh_table);
    const auto unique_id = bspline->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<RefSpaceType>(bspline);
  }
}



template <int dim>
void
ObjectsContainerXMLReader::
parse_identity_grid_function(const shared_ptr<XMLElement> xml_elem,
                             const bool parse_as_constant,
                             IdMap_ &id_map,
                             const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "IdentityGridFunction",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));

  using IdentityGridFunctionType = grid_functions::IdentityGridFunction<dim>;
  using GridFunctionType = GridFunction<dim, dim>;
  using GridType = Grid<dim>;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("IdentityGridFunction",
                                                       local_object_id, SafeSTLVector<int>(1, dim));

  const auto gr_tag = xml_elem->get_single_element("Grid");
  const auto local_gr_id = gr_tag->get_attribute<Index>("GetFromLocalObjectId");

  // Checking the grid with proper dimension and id exists.
  AssertThrow(id_map.find(local_gr_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<GridType> (id_map.at(local_gr_id)) :
               container->is_object_present<GridType> (id_map.at(local_gr_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching definition" +
                         "for " + Self_::get_type_id_string("Grid", local_gr_id,
                                                            SafeSTLVector<Index>(dim)) + "."));

  if (parse_as_constant)
  {
    const auto grid = container->get_const_object<GridType>(id_map.at(local_gr_id));
    const auto id_func = IdentityGridFunctionType::const_create(grid);
    const auto unique_id = id_func->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<GridFunctionType>(id_func);
  }
  else
  {
    const auto grid = container->get_object<GridType>(id_map.at(local_gr_id));
    const auto id_func = IdentityGridFunctionType::create(grid);
    const auto unique_id = id_func->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<GridFunctionType>(id_func);
  }
}



template <int dim, int space_dim>
void
ObjectsContainerXMLReader::
parse_constant_grid_function(const shared_ptr<XMLElement> xml_elem,
                             const bool parse_as_constant,
                             IdMap_ &id_map,
                             const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "ConstantGridFunction",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Spacedim") == space_dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Spacedim"), space_dim));

  using GridType = Grid<dim>;
  using ConstGridFunctionType = grid_functions::ConstantGridFunction<dim, space_dim>;
  using GridFunctionType = GridFunction<dim, space_dim>;
  using Values = typename GridFunctionType::Value;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("IgGridFunction",
  local_object_id, {{dim, space_dim}});

  // Gettting grid.
  const auto gr_tag = xml_elem->get_single_element("Grid");
  const auto local_gr_id = gr_tag->get_attribute<Index>("GetFromLocalObjectId");

  // Checking the grid with proper dimension and id exists.
  AssertThrow(id_map.find(local_gr_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<GridType> (id_map.at(local_gr_id)) :
               container->is_object_present<GridType> (id_map.at(local_gr_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching definition" +
                         "for " + Self_::get_type_id_string("Grid", local_gr_id,
                                                            SafeSTLVector<Index>(dim)) + "."));

  // Parsing values.
  const auto vals_tag = xml_elem->get_single_element("Values");

  const auto vals_vec = vals_tag->get_values_vector<Real>();
  AssertThrow(vals_vec.size() == Values::n_entries,
              ExcMessage("Parsing " + parsing_msg + ", the number of "
                         "components in Values XML does not match "
                         "with the number of components of the GridFunction."));
  SafeSTLArray<Real, Values::n_entries> vals_arr;
  std::copy(vals_vec.cbegin(), vals_vec.cend(), vals_arr.begin());
  Values values (vals_arr);



  if (parse_as_constant)
  {
    const auto grid = container->get_const_object<GridType>(id_map.at(local_gr_id));
    const auto cgf = ConstGridFunctionType::const_create(grid, values);
    const auto unique_id = cgf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<GridFunctionType>(cgf);
  }
  else
  {
    const auto grid = container->get_object<GridType>(id_map.at(local_gr_id));
    const auto cgf = ConstGridFunctionType::create(grid, values);
    const auto unique_id = cgf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<GridFunctionType>(cgf);
  }
}



template <int dim, int space_dim>
void
ObjectsContainerXMLReader::
parse_linear_grid_function(const shared_ptr<XMLElement> xml_elem,
                           const bool parse_as_constant,
                           IdMap_ &id_map,
                           const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "LinearGridFunction",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Spacedim") == space_dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Spacedim"), space_dim));

  using GridType = Grid<dim>;
  using LinearGridFunctionType = grid_functions::LinearGridFunction<dim, space_dim>;
  using GridFunctionType = GridFunction<dim, space_dim>;
  using Values = typename LinearGridFunctionType::Value;
  using Ders = typename LinearGridFunctionType::template Derivative<1>;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("IgGridFunction",
  local_object_id, {{dim, space_dim}});

  // Gettting grid.
  const auto gr_tag = xml_elem->get_single_element("Grid");
  const auto local_gr_id = gr_tag->get_attribute<Index>("GetFromLocalObjectId");

  // Checking the grid with proper dimension and id exists.
  AssertThrow(id_map.find(local_gr_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<GridType> (id_map.at(local_gr_id)) :
               container->is_object_present<GridType> (id_map.at(local_gr_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching definition" +
                         "for " + Self_::get_type_id_string("Grid", local_gr_id,
                                                            SafeSTLVector<Index>(dim)) + "."));

  // Parsing b.
  const auto b_tag = xml_elem->get_single_element("b");

  const auto b_vec = b_tag->get_values_vector<Real>();
  AssertThrow(b_vec.size() == Values::n_entries,
              ExcMessage("Parsing " + parsing_msg + ", the number of "
                         "components in Values XML does not match "
                         "with the number of components of the GridFunction."));
  SafeSTLArray<Real, Values::n_entries> b_arr;
  std::copy(b_vec.cbegin(), b_vec.cend(), b_arr.begin());
  Values b (b_arr);

  // Parsing A.
  const auto A_tag = xml_elem->get_single_element("A");

  const auto A_vec = A_tag->get_values_vector<Real>();
  AssertThrow(A_vec.size() == Ders::n_entries,
              ExcMessage("Parsing " + parsing_msg + ", the number of "
                         "components in Derivative<1> XML does not match "
                         "with the number of components of the Function."));
  SafeSTLArray<Real, Ders::n_entries> A_arr;
  std::copy(A_vec.cbegin(), A_vec.cend(), A_arr.begin());
  Ders A (A_arr);


  if (parse_as_constant)
  {
    const auto grid = container->get_const_object<GridType>(id_map.at(local_gr_id));
    const auto cgf = LinearGridFunctionType::const_create(grid, A, b);
    const auto unique_id = cgf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<GridFunctionType>(cgf);
  }
  else
  {
    const auto grid = container->get_object<GridType>(id_map.at(local_gr_id));
    const auto cgf = LinearGridFunctionType::create(grid, A, b);
    const auto unique_id = cgf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<GridFunctionType>(cgf);
  }
}



template <int dim, int space_dim>
void
ObjectsContainerXMLReader::
parse_ig_grid_function(const shared_ptr<XMLElement> xml_elem,
                       const bool parse_as_constant,
                       const bool &first_parsing,
                       IdMap_ &id_map,
                       const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "IgGridFunction",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Spacedim") == space_dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Spacedim"), space_dim));

  static const int rank = 1;
  using IgGridFunctionType = IgGridFunction<dim, space_dim>;
  using GridFunctionType = GridFunction<dim, space_dim>;
  using RefSpaceType     = ReferenceSpaceBasis<dim, space_dim, rank>;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");

  const bool already_parsed =
    id_map.find(local_object_id) != id_map.cend() &&
    (parse_as_constant ?
     container->is_const_object_present<GridFunctionType> (id_map.at(local_object_id)) :
     container->is_object_present<GridFunctionType> (id_map.at(local_object_id)));

  if (!first_parsing && already_parsed)
    return;

  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("IgGridFunction",
  local_object_id, {{dim, space_dim}});

  const auto rs_tag = xml_elem->get_single_element("ReferenceSpaceBasis");
  const auto local_rs_id = rs_tag->get_attribute<Index>("GetFromLocalObjectId");

  // Checking the reference space with proper dimension and id exists.
  const bool ref_space_parsed =
    id_map.find(local_rs_id) != id_map.cend() &&
    (parse_as_constant ?
     container->is_const_object_present<RefSpaceType> (id_map.at(local_rs_id)) :
     container->is_object_present<RefSpaceType> (id_map.at(local_rs_id)));

  if (!ref_space_parsed)
  {
    // If does not exist there are two possibiities:
    if (first_parsing)
      // It is the first parsing time for the ig grid function.
      // It is assumed that the reference space is a NURBS that
      // has not been parsed yet.
      return;
    else
      // It is not the first parsing time for the ig grid function.
      // There is an error in the input file
      AssertThrow(false,
                  ExcMessage("Parsing " + parsing_msg + " not matching "
                             "definition for " +
                             Self_::get_type_id_string("ReferenceSpaceBasis", local_rs_id,
    {{dim, space_dim, rank}}) + "."));

  }

  SharedPtrConstnessHandler<RefSpaceType> rs;
  if (parse_as_constant)
    rs = SharedPtrConstnessHandler<RefSpaceType>
         (container->get_const_object<RefSpaceType>(id_map.at(local_rs_id)));
  else
    rs = SharedPtrConstnessHandler<RefSpaceType>
         (container->get_object<RefSpaceType>(id_map.at(local_rs_id)));

  const string dofs_property = parse_dofs_property(xml_elem);
  const auto dof_distribution = rs->get_spline_space()->get_dof_distribution();
  AssertThrow(dof_distribution->is_property_defined(dofs_property),
              ExcMessage("Parsing " + parsing_msg + " dofs property \"" +
                         dofs_property + "\" not defined for " +
                         Self_::get_type_id_string("ReferenceSpaceBasis", local_rs_id,
  {{dim, space_dim, rank}}) + "."));

  const auto &global_dofs = dof_distribution->get_global_dofs(dofs_property);
  const auto ig_coefs = parse_ig_coefficients(xml_elem, parsing_msg, global_dofs);

  AssertThrow(rs->get_num_basis() == ig_coefs->size(),
              ExcMessage("Parsing " + parsing_msg + " the cardinality "
                         "of the ReferenceSpaceBasis (" +
                         to_string(rs->get_num_basis()) + ") is "
                         "different to the dimension of the "
                         "IgCoefficients (" + to_string(ig_coefs->size())
                         + ")."));

  if (parse_as_constant)
  {
    const auto igf = IgGridFunctionType::const_create
                     (rs.get_ptr_const_data(), *ig_coefs, dofs_property);
    const auto unique_id = igf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<GridFunctionType>(igf);
  }
  else
  {
    const auto igf = IgGridFunctionType::create(rs.get_ptr_data(), *ig_coefs, dofs_property);
    const auto unique_id = igf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<GridFunctionType>(igf);
  }
}



template <int dim, int range, int rank>
void
ObjectsContainerXMLReader::
parse_nurbs(const shared_ptr<XMLElement> xml_elem,
            const bool parse_as_constant,
            IdMap_ &id_map,
            const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "NURBS",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Range") == range,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
  Assert(xml_elem->get_attribute<int>("Rank") == rank,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));

  using BSplineType = BSpline<dim, range, rank>;
  using NURBSType   = NURBS<dim, range, rank>;
  using WeightIgFunctionType = IgGridFunction<dim, 1>;
  using WeightFunctionType = GridFunction<dim, 1>;
  using RefSpaceType = ReferenceSpaceBasis<dim, range, rank>;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("NURBS",
  local_object_id, {{dim, range, rank}});

  // Parsing BSpline
  const auto bs_tag = xml_elem->get_single_element("BSpline");
  const auto local_bs_id = bs_tag->get_attribute<Index>("GetFromLocalObjectId");
  // Checking the bspline with the proper dimensions and id exists.
  AssertThrow(id_map.find(local_bs_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<RefSpaceType> (id_map.at(local_bs_id)) :
               container->is_object_present<RefSpaceType> (id_map.at(local_bs_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("BSpline", local_bs_id,
  {{dim, range, rank}}) + "."));

  SharedPtrConstnessHandler<RefSpaceType> rs;
  if (parse_as_constant)
    rs = SharedPtrConstnessHandler<RefSpaceType>
         (container->get_const_object<RefSpaceType>(id_map.at(local_bs_id)));
  else
    rs = SharedPtrConstnessHandler<RefSpaceType>
         (container->get_object<RefSpaceType>(id_map.at(local_bs_id)));

  // Checking that the reference space is a BSpline.
  AssertThrow(rs->is_bspline(),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("BSplineSpaceBasis", local_bs_id,
  {{dim, range, rank}}) + ". It is a "
  "ReferenceBasisSpace, but not a BSpline."
  "BSpline space basis must be used for "
  "defining IgGridFunctions used as weight "
  "functions for NURBS space basis."));

  // Parsing Weight function
  const auto wg_tag = xml_elem->get_single_element("WeightFunction");
  const auto local_wg_id = wg_tag->get_attribute<Index>("GetFromLocalObjectId");
  // Checking that the grid function with the proper dimensions exists.
  AssertThrow(id_map.find(local_wg_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<WeightFunctionType> (id_map.at(local_wg_id)) :
               container->is_object_present<WeightFunctionType> (id_map.at(local_wg_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("IgGridFunction", local_wg_id,
  {{dim, 1}}) +  ". The error could be caused "
  "by the fact the object does not correspond "
  "to IgGridFunction built upon a BSpline. "
  "Currently igatools only allows to build "
  "weight functions based on scalar BSpline "
  "space basis."));

  SharedPtrConstnessHandler<WeightIgFunctionType> wf;
  if (parse_as_constant)
  {
    const auto gf = container->get_const_object<WeightFunctionType>(id_map.at(local_wg_id));
    wf = SharedPtrConstnessHandler<WeightIgFunctionType>(
           dynamic_pointer_cast<const WeightIgFunctionType>(gf));
  }
  else
  {
    const auto gf = container->get_object<WeightFunctionType>(id_map.at(local_wg_id));
    wf = SharedPtrConstnessHandler<WeightIgFunctionType>(
           dynamic_pointer_cast<WeightIgFunctionType>(gf));
  }

  // Checking that the grid function is a ig grid function.
  AssertThrow((parse_as_constant && wf.get_ptr_const_data()) ||
              (!parse_as_constant && wf.get_ptr_data()),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("IgGridFunction", local_wg_id,
  {{dim, 1}}) + ". It is a"
  "GridFunction, but not a IgGridFunction"));

  // Checking that the weight function was defined upon a BSpline.
  const auto w_func_basis = wf->get_basis();
  Assert(wf->get_basis()->is_bspline(),
         ExcMessage("Parsing " + parsing_msg + ", " +
                    Self_::get_type_id_string("IgGridFunction", local_wg_id,
  {{dim, 1}}) + ". It is based on a NURBS space basis"
  ", but must be based on BSpline space basis."));

  // Checking that the grids of the bspline and weight function match.
  AssertThrow(rs->get_grid() == wf->get_grid(),
              ExcMessage("Parsing " + parsing_msg + ", mismatching "
                         "grids for the BSpline and WeightFunction."));

  // Checking that the number of basis functions and weights match.
  const auto &n_basis_table = rs->get_spline_space()->get_dof_distribution()->get_num_dofs_table();
  int comp_id = 0;
  for (const auto &n_basis_comp : n_basis_table)
  {
    AssertThrow(n_basis_comp ==
                w_func_basis->get_spline_space()->get_dof_distribution()->get_num_dofs_table()[0],
                ExcMessage("Parsing " + parsing_msg + ", mismatching number of basis functions and weight "
                           "coefficients for scalar component " + to_string(comp_id++)));
  }

  if (parse_as_constant)
  {
    const auto bs = dynamic_pointer_cast<const BSplineType>(rs.get_ptr_const_data());
    const auto nurbs = NURBSType::const_create(bs, wf.get_ptr_const_data());
    Assert(nurbs != nullptr, ExcNullPtr());
    const auto unique_id = nurbs->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<RefSpaceType>(nurbs);
  }
  else
  {
    const auto bs = dynamic_pointer_cast<BSplineType>(rs.get_ptr_data());
    const auto nurbs = NURBSType::create(bs, wf.get_ptr_data());
    Assert(nurbs != nullptr, ExcNullPtr());
    const auto unique_id = nurbs->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<RefSpaceType>(nurbs);
  }
}



template <int dim, int codim>
void
ObjectsContainerXMLReader::
parse_domain(const shared_ptr<XMLElement> xml_elem,
             const bool parse_as_constant,
             IdMap_ &id_map,
             const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "Domain",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Codim") == codim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Codim"), codim));

  static const int space_dim = dim + codim;
  using DomainType = Domain<dim, codim>;
  using GridFuncType = GridFunction<dim, space_dim>;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("Domain",
  local_object_id, {{dim, codim}});

  const auto gf_tag = xml_elem->get_single_element("GridFunction");
  const auto local_gf_id = gf_tag->get_attribute<Index>("GetFromLocalObjectId");
  // Checking if the grid function exists.
  AssertThrow(id_map.find(local_gf_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<GridFuncType> (id_map.at(local_gf_id)) :
               container->is_object_present<GridFuncType> (id_map.at(local_gf_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("GridFunction", local_gf_id,
  {{dim, space_dim}}) + "."));

  const auto name = parse_name(xml_elem);

  if (parse_as_constant)
  {
    const auto gf = container->get_const_object<GridFuncType>(id_map.at(local_gf_id));
    const auto domain = DomainType::const_create(gf, name);
    const auto unique_id = domain->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<DomainType>(domain);
  }
  else
  {
    const auto gf = container->get_object<GridFuncType>(id_map.at(local_gf_id));
    const auto domain = DomainType::create(gf, name);
    const auto unique_id = domain->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<DomainType>(domain);
  }
}



template <int dim, int codim, int range, int rank>
void
ObjectsContainerXMLReader::
parse_ig_function(const shared_ptr<XMLElement> xml_elem,
                  const bool parse_as_constant,
                  IdMap_ &id_map,
                  const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "IgFunction",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Range") == range,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
  Assert(xml_elem->get_attribute<int>("Rank") == rank,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));
  Assert(xml_elem->get_attribute<int>("Codim") == codim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Codim"), codim));

  using PhysSpaceType = PhysicalSpaceBasis<dim, range, rank, codim>;
  using IgFunctionType = IgFunction<dim, codim, range, rank>;
  using FunctionType = Function<dim, codim, range, rank>;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("IgFunction",
  local_object_id, {{dim, codim, range, rank}});

  const auto ps_tag = xml_elem->get_single_element("PhysicalSpaceBasis");
  const auto local_ps_id = ps_tag->get_attribute<Index>("GetFromLocalObjectId");

  // Checking if the physical space exists.
  AssertThrow(id_map.find(local_ps_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<PhysSpaceType> (id_map.at(local_ps_id)) :
               container->is_object_present<PhysSpaceType> (id_map.at(local_ps_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("PhysicalSpaceBasis", local_ps_id,
  {{dim, range, rank, codim}}) + "."));

  SharedPtrConstnessHandler<PhysSpaceType> ps;
  if (parse_as_constant)
    ps = SharedPtrConstnessHandler<PhysSpaceType>
         (container->get_const_object<PhysSpaceType>(id_map.at(local_ps_id)));
  else
    ps = SharedPtrConstnessHandler<PhysSpaceType>
         (container->get_object<PhysSpaceType>(id_map.at(local_ps_id)));

  const auto name = parse_name(xml_elem);

  const string dofs_property = parse_dofs_property(xml_elem);
  const auto dof_distribution = ps->get_spline_space()->get_dof_distribution();
  AssertThrow(dof_distribution->is_property_defined(dofs_property),
              ExcMessage("Parsing " + parsing_msg + " dofs property \"" +
                         dofs_property + "\" not defined for " +
                         Self_::get_type_id_string("PhysicalSpaceBasis", local_ps_id,
  {{dim, range, rank, codim}}) + "."));

  const auto &global_dofs = dof_distribution->get_global_dofs(dofs_property);
  const auto ig_coefs = parse_ig_coefficients(xml_elem, parsing_msg, global_dofs);

  AssertThrow(ps->get_num_basis() == ig_coefs->size(),
              ExcMessage("Parsing " + parsing_msg + " the cardinality "
                         "of the PhysicalSpaceBasis (" +
                         to_string(ps->get_num_basis()) + ") is "
                         "different to the dimension of the "
                         "IgCoefficients (" + to_string(ig_coefs->size())
                         + ")."));

  if (parse_as_constant)
  {
    const auto igf = IgFunctionType::const_create(ps.get_ptr_const_data(), *ig_coefs, dofs_property, name);
    const auto unique_id = igf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<FunctionType>(igf);
  }
  else
  {
    const auto igf = IgFunctionType::create(ps.get_ptr_data(), *ig_coefs, dofs_property, name);
    const auto unique_id = igf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<FunctionType>(igf);
  }
}



template <int dim, int codim, int range, int rank>
void
ObjectsContainerXMLReader::
parse_constant_function(const shared_ptr<XMLElement> xml_elem,
                        const bool parse_as_constant,
                        IdMap_ &id_map,
                        const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "ConstantFunction",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Range") == range,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
  Assert(xml_elem->get_attribute<int>("Rank") == rank,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));
  Assert(xml_elem->get_attribute<int>("Codim") == codim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Codim"), codim));

  using DomainType = Domain<dim, codim>;
  using FunctionType = Function<dim, codim, range, rank>;
  using ConstFunctionType = functions::ConstantFunction<dim, codim, range, rank>;
  using Values = typename FunctionType::Value;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("ConstantFunction",
  local_object_id, {{dim, codim, range, rank}});

  const auto dm_tag = xml_elem->get_single_element("Domain");
  const auto local_dm_id = dm_tag->get_attribute<Index>("GetFromLocalObjectId");

  AssertThrow(id_map.find(local_dm_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<DomainType> (id_map.at(local_dm_id)) :
               container->is_object_present<DomainType> (id_map.at(local_dm_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("Domain", local_dm_id,
  {{dim, codim}}) + "."));

  SharedPtrConstnessHandler<DomainType> dm;
  if (parse_as_constant)
    dm = SharedPtrConstnessHandler<DomainType>
         (container->get_const_object<DomainType>(id_map.at(local_dm_id)));
  else
    dm = SharedPtrConstnessHandler<DomainType>
         (container->get_object<DomainType>(id_map.at(local_dm_id)));

  // Parsing values.
  const auto vals_tag = xml_elem->get_single_element("Values");

  const auto vals_vec = vals_tag->get_values_vector<Real>();
  AssertThrow(vals_vec.size() == Values::n_entries,
              ExcMessage("Parsing " + parsing_msg + ", the number of "
                         "components in Values XML does not match "
                         "with the number of components of the GridFunction."));
  SafeSTLArray<Real, Values::n_entries> vals_arr;
  std::copy(vals_vec.cbegin(), vals_vec.cend(), vals_arr.begin());
  Values values (vals_arr);

  const auto name = parse_name(xml_elem);

  if (parse_as_constant)
  {
    const auto cf = ConstFunctionType::const_create(dm.get_ptr_const_data(), values, name);
    const auto unique_id = cf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<FunctionType>(cf);
  }
  else
  {
    const auto cf = ConstFunctionType::create(dm.get_ptr_data(), values, name);
    const auto unique_id = cf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<FunctionType>(cf);
  }
}



template <int dim, int codim, int range>
void
ObjectsContainerXMLReader::
parse_linear_function(const shared_ptr<XMLElement> xml_elem,
                      const bool parse_as_constant,
                      IdMap_ &id_map,
                      const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "LinearFunction",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Range") == range,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
  Assert(xml_elem->get_attribute<int>("Codim") == codim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Codim"), codim));

  static const int rank = 1;
  using DomainType = Domain<dim, codim>;
  using FunctionType = Function<dim, codim, range, rank>;
  using LinearFunctionType = functions::LinearFunction<dim, codim, range>;
  using Values = typename FunctionType::Value;
  using Ders = typename FunctionType::template Derivative<1>;

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const string parsing_msg = Self_::get_type_id_string("LinearFunction",
  local_object_id, {{dim, codim, range, rank}});

  const auto dm_tag = xml_elem->get_single_element("Domain");
  const auto local_dm_id = dm_tag->get_attribute<Index>("GetFromLocalObjectId");

  AssertThrow(id_map.find(local_dm_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<DomainType> (id_map.at(local_dm_id)) :
               container->is_object_present<DomainType> (id_map.at(local_dm_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("Domain", local_dm_id,
  {{dim, codim}}) + "."));

  SharedPtrConstnessHandler<DomainType> dm;
  if (parse_as_constant)
    dm = SharedPtrConstnessHandler<DomainType>
         (container->get_const_object<DomainType>(id_map.at(local_dm_id)));
  else
    dm = SharedPtrConstnessHandler<DomainType>
         (container->get_object<DomainType>(id_map.at(local_dm_id)));

  // Parsing b.
  const auto b_tag = xml_elem->get_single_element("b");

  const auto b_vec = b_tag->get_values_vector<Real>();
  AssertThrow(b_vec.size() == Values::n_entries,
              ExcMessage("Parsing " + parsing_msg + ", the number of "
                         "components in Values XML does not match "
                         "with the number of components of the Function."));
  SafeSTLArray<Real, Values::n_entries> b_arr;
  std::copy(b_vec.cbegin(), b_vec.cend(), b_arr.begin());
  Values b (b_arr);

  // Parsing A.
  const auto A_tag = xml_elem->get_single_element("A");

  const auto A_vec = A_tag->get_values_vector<Real>();
  AssertThrow(A_vec.size() == Ders::n_entries,
              ExcMessage("Parsing " + parsing_msg + ", the number of "
                         "components in Derivative<1> XML does not match "
                         "with the number of components of the Function."));
  SafeSTLArray<Real, Ders::n_entries> A_arr;
  std::copy(A_vec.cbegin(), A_vec.cend(), A_arr.begin());
  Ders A (A_arr);

  const auto name = parse_name(xml_elem);

  if (parse_as_constant)
  {
    const auto cf = LinearFunctionType::const_create(dm.get_ptr_const_data(), A, b, name);
    const auto unique_id = cf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<FunctionType>(cf);
  }
  else
  {
    const auto cf = LinearFunctionType::create(dm.get_ptr_data(), A, b, name);
    const auto unique_id = cf->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<FunctionType>(cf);
  }
}



template <int dim, int codim, int range, int rank>
void
ObjectsContainerXMLReader::
parse_phys_space(const shared_ptr<XMLElement> xml_elem,
                 const bool parse_as_constant,
                 IdMap_ &id_map,
                 const shared_ptr<ObjectsContainer> container)
{
  Assert(xml_elem->get_name() == "PhysicalSpaceBasis",
         ExcMessage("Invalid XML tag."));

  Assert(xml_elem->get_attribute<int>("Dim") == dim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
  Assert(xml_elem->get_attribute<int>("Range") == range,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
  Assert(xml_elem->get_attribute<int>("Rank") == rank,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));
  Assert(xml_elem->get_attribute<int>("Codim") == codim,
         ExcDimensionMismatch(xml_elem->get_attribute<int>("Codim"), codim));

  using PhysSpaceType = PhysicalSpaceBasis<dim, range, rank, codim>;
  using RefSpaceType = ReferenceSpaceBasis<dim, range, rank>;
  using DomainType = Domain<dim, codim>;

  Transformation transf = Transformation::h_grad;
  if (xml_elem->has_element("Transformation"))
  {
    const auto trans_str = xml_elem->get_single_element("Transformation")->get_value<string>();
    if (trans_str == "h_grad")
      transf = Transformation::h_grad;
    else if (trans_str == "h_div")
      transf = Transformation::h_div;
    else if (trans_str == "h_curl")
      transf = Transformation::h_curl;
    else if (trans_str == "l_2")
      transf = Transformation::l_2;
    // If is not one of those types, the XML schema should have
    // thrown the error before.
  }

  const auto local_object_id = xml_elem->get_attribute<Index>("LocalObjectId");
  Assert(id_map.find(local_object_id) == id_map.cend(), ExcMessage("Repeated object id."));

  const auto parsing_msg = Self_::get_type_id_string("PhysicalSpaceBasis",
  local_object_id, {{dim, range, rank}});

  const auto rs_tag = xml_elem->get_single_element("ReferenceSpaceBasis");
  const auto local_rs_id = rs_tag->get_attribute<Index>("GetFromLocalObjectId");

  AssertThrow(id_map.find(local_rs_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<RefSpaceType> (id_map.at(local_rs_id)) :
               container->is_object_present<RefSpaceType> (id_map.at(local_rs_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("ReferenceSpaceBasis", local_rs_id,
  {{dim, range, rank}}) + "."));

  SharedPtrConstnessHandler<RefSpaceType> rs;
  if (parse_as_constant)
    rs = SharedPtrConstnessHandler<RefSpaceType>
         (container->get_const_object<RefSpaceType>(id_map.at(local_rs_id)));
  else
    rs = SharedPtrConstnessHandler<RefSpaceType>
         (container->get_object<RefSpaceType>(id_map.at(local_rs_id)));

  const auto dm_tag = xml_elem->get_single_element("Domain");
  const auto local_dm_id = dm_tag->get_attribute<Index>("GetFromLocalObjectId");

  AssertThrow(id_map.find(local_dm_id) != id_map.cend() &&
              (parse_as_constant ?
               container->is_const_object_present<DomainType> (id_map.at(local_dm_id)) :
               container->is_object_present<DomainType> (id_map.at(local_dm_id))),
              ExcMessage("Parsing " + parsing_msg + " not matching "
                         "definition for " +
                         Self_::get_type_id_string("Domain", local_dm_id,
  {{dim, codim}}) + "."));

  SharedPtrConstnessHandler<DomainType> dm;
  if (parse_as_constant)
    dm = SharedPtrConstnessHandler<DomainType>
         (container->get_const_object<DomainType>(id_map.at(local_dm_id)));
  else
    dm = SharedPtrConstnessHandler<DomainType>
         (container->get_object<DomainType>(id_map.at(local_dm_id)));

  // Checking that the reference space and domain grids match.
  AssertThrow(rs->get_grid() == dm->get_grid_function()->get_grid(),
              ExcMessage("Parsing " + parsing_msg + ", mismatching "
                         "grids for the ReferenceSpaceBasis and Domain."));

  if (parse_as_constant)
  {
    const auto ps = PhysSpaceType::const_create(rs.get_ptr_const_data(), dm.get_ptr_const_data(), transf);
    const auto unique_id = ps->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_const_object<PhysSpaceType>(ps);
  }
  else
  {
    const auto ps = PhysSpaceType::create(rs.get_ptr_data(), dm.get_ptr_data(), transf);
    const auto unique_id = ps->get_object_id();
    id_map[local_object_id] = unique_id;

    container->insert_object<PhysSpaceType>(ps);
  }
}



string
ObjectsContainerXMLReader::
parse_name(const shared_ptr<XMLElement> xml_elem)
{
  if (xml_elem->has_element("Name"))
    return xml_elem->get_single_element("Name")->get_value<string>();
  else
    return string("");
}



string
ObjectsContainerXMLReader::
parse_dofs_property(const shared_ptr<XMLElement> xml_elem)
{
  if (xml_elem->has_element("DofsProperty"))
    return xml_elem->get_single_element("DofsProperty")->get_value<string>();
  else
    return DofProperties::active;
}



string
ObjectsContainerXMLReader::
get_type_id_string(const string &object_type,
                   const Index &local_object_id,
                   const SafeSTLVector<int> &dims)
{

  string dims_str = object_type + "<";
  for (const auto &d : dims)
    dims_str += to_string(d) + ", ";
  return dims_str.substr(0, dims_str.size() - 2) + ">" +
         " (LocalObjectId=" + to_string(local_object_id) + ")";
}



shared_ptr<IgCoefficients>
ObjectsContainerXMLReader::
parse_ig_coefficients(const shared_ptr<XMLElement> xml_elem,
                      const string &parsing_msg,
                      const set<Index> &space_global_dofs)
{
  Assert(xml_elem->has_element("IgCoefficients"),
         ExcMessage("IgCoefficients XML element not present."));

  const auto ig_elem = xml_elem->get_single_element("IgCoefficients");
  const auto size = ig_elem->get_attribute<Index>("Size");

  // Gettint the indices.
  const auto ig_ind_elem = ig_elem->get_single_element("Indices");
  const auto ig_coefs_ind_vec = ig_ind_elem->get_values_vector<Index>();
  AssertThrow(size == ig_coefs_ind_vec.size(),
              ExcMessage("Parsing IgCoefficients indices for " + parsing_msg +
                         ", Size=" + to_string(size) + " do not match "
                         "with the vector size."));

  // Gettint the values.
  const auto ig_val_elem = ig_elem->get_single_element("Values");
  const auto ig_coefs_val_vec = ig_val_elem->get_values_vector<Real>();
  AssertThrow(size == ig_coefs_val_vec.size(),
              ExcMessage("Parsing IgCoefficients values for " + parsing_msg +
                         ", Size=" + to_string(size) + " do not match "
                         "with the vector size."));

  set<Index> indices;
  for (const auto &i : ig_coefs_ind_vec)
    indices.insert(i);
  copy(ig_coefs_ind_vec.cbegin(), ig_coefs_ind_vec.cend(),
       inserter(indices, indices.begin()));

  // Checking that there are not repeated indices.
  AssertThrow(size == indices.size(),
              ExcMessage("Parsing IgCoefficients for " + parsing_msg +
                         ", not valid indices vector parsed. Repeated "
                         "indices may found."));

  const auto ig_coefs = shared_ptr<IgCoefficients>(new IgCoefficients(indices));

  // Checking if the parsed indices match with space_global_dofs and
  // filling the values of the ig coefficients vector.
  const auto end_dofs = space_global_dofs.cend();
  auto &igc = *ig_coefs;
  auto ind_it = ig_coefs_ind_vec.cbegin();
  for (const auto &val : ig_coefs_val_vec)
  {
    AssertThrow(space_global_dofs.find(*ind_it) != end_dofs,
                ExcMessage("Parsing IgCoefficients for " + parsing_msg +
                           ", " + to_string(*ind_it) + " is not a valid index."));
    igc[*ind_it++] = val;
  }

  return ig_coefs;
}

IGA_NAMESPACE_CLOSE

#endif // XML_IO
