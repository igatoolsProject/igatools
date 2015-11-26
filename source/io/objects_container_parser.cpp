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

#include <igatools/io/objects_container_parser.h>

#ifdef XML_IO

#include <igatools/io/xml_file_parser.h>
#include <igatools/io/xml_element.h>
#include <igatools/base/objects_container.h>
#include <igatools/base/instantiated_types.h>

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
ObjectsContainerParser::
parse(const string &file_path)
{
    const auto parser = XMLFileParser::create();
    const string schema_file = "include/igatools/io/objects_container_XML_schema.xsd";
    const auto xml_elem = parser->parse(file_path, schema_file);
    const auto container = ObjectsContainer::create();

    // Checking for repeated iga object ids.
    SafeSTLSet<Index> object_ids;
    const auto kk = xml_elem->get_children_elements();
    for (const auto &el : xml_elem->get_children_elements())
    {
        const Index obj_id = el->get_attribute<Index>("IgaObjectId");
        AssertThrow (object_ids.find(obj_id) == object_ids.cend(),
                ExcMessage("IgaObjectId " + to_string(obj_id) + " is "
                           "used more than once."));
        object_ids.insert(obj_id);
    }

    Self_::parse_grids (xml_elem, container);
    Self_::parse_spline_spaces (xml_elem, container);
    Self_::parse_bsplines (xml_elem, container);
    Self_::parse_grid_functions_and_nurbs (xml_elem, container);
    Self_::parse_domains (xml_elem, container);
    Self_::parse_phys_spaces (xml_elem, container);
    Self_::parse_functions (xml_elem, container);

    return container;
}



void
ObjectsContainerParser::
parse_grids(const shared_ptr<XMLElement> xml_elem,
            const shared_ptr<ObjectsContainer> container)
{
    const auto grid_elems = xml_elem->get_children_elements("Grid");
    for (const auto &ge : grid_elems)
    {
        const int grid_dim = ge->get_attribute<int>("Dim");

        using GridPtrs = typename InstantiatedTypes::GridPtrs;
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
                parse_grid<dim>(ge, container);
            }
        });

        // Grid dimension not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("Grid",
                     ge->get_attribute<Index>("IgaObjectId"),
                     SafeSTLVector<int>(1, grid_dim))
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_spline_spaces(const shared_ptr<XMLElement> xml_elem,
                    const shared_ptr<ObjectsContainer> container)
{
    for (const auto &ssp : xml_elem->get_children_elements("SplineSpace"))
    {
        const int ssp_dim = ssp->get_attribute<int>("Dim");
        const int ssp_range = ssp->get_attribute<int>("Range");
        const int ssp_rank = ssp->get_attribute<int>("Rank");

        using SplineSpacePtrs = typename InstantiatedTypes::SplineSpacePtrs;
        SplineSpacePtrs valid_ssp_ptr_types;

        bool found = false;
        boost::fusion::for_each(valid_ssp_ptr_types, [&](const auto &ssp_ptr_type)
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
                parse_spline_space<dim, range, rank>(ssp, container);
            }
        });

        // SplineSpace dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("SplineSpace",
                     ssp->get_attribute<Index>("IgaObjectId"),
                     {{ssp_dim, ssp_range, ssp_rank}})
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_bsplines(const shared_ptr<XMLElement> xml_elem,
               const shared_ptr<ObjectsContainer> container)
{
    for (const auto &bs : xml_elem->get_children_elements("BSpline"))
    {
        const int bs_dim = bs->get_attribute<int>("Dim");
        const int bs_range = bs->get_attribute<int>("Range");
        const int bs_rank = bs->get_attribute<int>("Rank");

        using BSplinePtrs = typename InstantiatedTypes::SplineSpacePtrs;
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
                parse_bspline<dim, range, rank>(bs, container);
            }
        });

        // BSpline dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("BSpline",
                     bs->get_attribute<Index>("IgaObjectId"),
                     {{bs_dim, bs_range, bs_rank}})
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_nurbs(const shared_ptr<XMLElement> xml_elem,
            const shared_ptr<ObjectsContainer> container)
{
    for (const auto &nr : xml_elem->get_children_elements("NURBS"))
    {
        const int nr_dim = nr->get_attribute<int>("Dim");
        const int nr_range = nr->get_attribute<int>("Range");
        const int nr_rank = nr->get_attribute<int>("Rank");

        using NURBSPtrs = typename InstantiatedTypes::SplineSpacePtrs;
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
                parse_nurbs<dim, range, rank>(nr, container);
            }
        });

        // NURBS dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("NURBS",
                     nr->get_attribute<Index>("IgaObjectId"),
                     {{nr_dim, nr_range, nr_rank}})
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_identity_grid_functions(const shared_ptr<XMLElement> xml_elem,
                              const shared_ptr<ObjectsContainer> container)
{
    for (const auto &id : xml_elem->get_children_elements("IdentityGridFunction"))
    {
        const int id_dim = id->get_attribute<int>("Dim");

        using GridPtrs = typename InstantiatedTypes::GridPtrs;
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
                parse_identity_grid_function<dim>(id, container);
            }
        });

        // IdentityGridFunction dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("IdentityGridFunction",
                     id->get_attribute<Index>("IgaObjectId"),
                     SafeSTLVector<int>(1, id_dim))
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_constant_grid_functions(const shared_ptr<XMLElement> xml_elem,
                              const shared_ptr<ObjectsContainer> container)
{
    for (const auto &cgf : xml_elem->get_children_elements("ConstantGridFunction"))
    {
        const int cgf_dim = cgf->get_attribute<int>("Dim");
        const int cgf_space_dim = cgf->get_attribute<int>("Spacedim"); // This is going to change.

        using GridFunctionPtrs = typename InstantiatedTypes::GridFunctionPtrs;
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
                parse_constant_grid_function<dim, space_dim>(cgf, container);
            }
        });

        // ConstantGridFunction dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("ConstantGridFunction",
                     cgf->get_attribute<Index>("IgaObjectId"),
                     {{cgf_dim, cgf_space_dim}})
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_ig_grid_functions(const shared_ptr<XMLElement> xml_elem,
                        const bool &first_parsing,
                        const shared_ptr<ObjectsContainer> container)
{
    for (const auto &gf : xml_elem->get_children_elements("IgGridFunction"))
    {
        const int gf_dim = gf->get_attribute<int>("Dim");
        const int gf_space_dim = gf->get_attribute<int>("Spacedim"); // This is going to change.

        using GridFunctionPtrs = typename InstantiatedTypes::GridFunctionPtrs;
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
                parse_ig_grid_function<dim, space_dim>(gf, first_parsing, container);
            }
        });

        // NURBS dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("IgGridFunction",
                     gf->get_attribute<Index>("IgaObjectId"),
                     {{gf_dim, gf_space_dim}})
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_grid_functions_and_nurbs(const shared_ptr<XMLElement> xml_elem,
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
    parse_identity_grid_functions(xml_elem, container);

    // Parsing constant grid functions.
    parse_constant_grid_functions(xml_elem, container);

    // Parsing ig grid functions built upon a BSpline.
    bool first_parsing = true;
    parse_ig_grid_functions(xml_elem, first_parsing, container);

    // Parsing NURBS.
    Self_::parse_nurbs(xml_elem, container);

    // Parsing the remaining ig grid functions.
    first_parsing = false;
    parse_ig_grid_functions(xml_elem, first_parsing, container);
}



void
ObjectsContainerParser::
parse_domains(const shared_ptr<XMLElement> xml_elem,
              const shared_ptr<ObjectsContainer> container)
{
    for (const auto &dm : xml_elem->get_children_elements("Domain"))
    {
        const int dm_dim = dm->get_attribute<int>("Dim");
        const int dm_codim = dm->get_attribute<int>("Codim");

        using DomainPtrs = typename InstantiatedTypes::DomainPtrs;
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
                parse_domain<dim, codim>(dm, container);
            }
        });

        // Domains dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("Domain",
                     dm->get_attribute<Index>("IgaObjectId"),
                     {{dm_dim, dm_codim}})
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_phys_spaces(const shared_ptr<XMLElement> xml_elem,
                  const shared_ptr<ObjectsContainer> container)
{
    for (const auto &ps : xml_elem->get_children_elements("PhysicalSpaceBasis"))
    {
        const int ps_dim = ps->get_attribute<int>("Dim");
        const int ps_codim = ps->get_attribute<int>("Codim");
        const int ps_range = ps->get_attribute<int>("Range");
        const int ps_rank = ps->get_attribute<int>("Rank");

        using PhysSpacePtrs = typename InstantiatedTypes::PhysSpacePtrs;
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
                parse_phys_space<dim, codim, range, rank>(ps, container);
            }
        });

        // PhysicalSpaceBasis dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("PhysicalSpaceBasis",
                     ps->get_attribute<Index>("IgaObjectId"),
                     {{ps_dim, ps_range, ps_rank, ps_codim}})
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_functions(const shared_ptr<XMLElement> xml_elem,
                const shared_ptr<ObjectsContainer> container)
{
    parse_ig_functions (xml_elem, container);
    parse_constant_functions (xml_elem, container);
}



void
ObjectsContainerParser::
parse_ig_functions(const shared_ptr<XMLElement> xml_elem,
                   const shared_ptr<ObjectsContainer> container)
{
    for (const auto &fn : xml_elem->get_children_elements("IgFunction"))
    {
        const int fn_dim = fn->get_attribute<int>("Dim");
        const int fn_codim = fn->get_attribute<int>("Codim");
        const int fn_range = fn->get_attribute<int>("Range");
        const int fn_rank = fn->get_attribute<int>("Rank");

        using FunctionPtrs = typename InstantiatedTypes::FunctionPtrs;
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
                parse_ig_function<dim, codim, range, rank>(fn, container);
            }
        });

        // Function dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("IgFunction",
                     fn->get_attribute<Index>("IgaObjectId"),
                     {{fn_dim, fn_codim, fn_range, fn_rank}})
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



void
ObjectsContainerParser::
parse_constant_functions(const shared_ptr<XMLElement> xml_elem,
                   const shared_ptr<ObjectsContainer> container)
{
    for (const auto &fn : xml_elem->get_children_elements("ConstantFunction"))
    {
        const int fn_dim = fn->get_attribute<int>("Dim");
        const int fn_codim = fn->get_attribute<int>("Codim");
        const int fn_range = fn->get_attribute<int>("Range");
        const int fn_rank = fn->get_attribute<int>("Rank");

        using FunctionPtrs = typename InstantiatedTypes::FunctionPtrs;
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
                parse_constant_function<dim, codim, range, rank>(fn, container);
            }
        });

        // Function dimensions not found
        AssertThrow (found,
          ExcMessage(Self_::get_type_id_string("ConstantFunction",
                     fn->get_attribute<Index>("IgaObjectId"),
                     {{fn_dim, fn_codim, fn_range, fn_rank}})
                     + " is not a valid type. Possibly the type was not "
                     "instantiated for the specified dimensions."));
    }
}



template <int dim>
void
ObjectsContainerParser::
parse_grid(const shared_ptr<XMLElement> xml_elem,
           const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "Grid", ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));

    using GridType = Grid<dim>;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    const string parsing_msg = Self_::get_type_id_string("Grid", object_id,
                               SafeSTLVector<int>(1, dim));

    const auto knots_children = xml_elem->get_children_elements("Knots");
    // Checking the number of knot vector match with the dimension.
    AssertThrow (dim == knots_children.size(),
                 ExcMessage("Parsing " + parsing_msg + ", the number of "
                            "knots vectors is not valid."));

    SafeSTLArray<SafeSTLVector<Real>, dim> knots;

    SafeSTLSet<Index> parsed_dirs;

    for (const auto &ke : knots_children)
    {
        const auto dir = ke->get_attribute<Index>("Direction");
        // Checking the <= 0 direction < dim
        AssertThrow (dir >= 0 && dir < dim,
                     ExcMessage("Parsing knot vectors for " + parsing_msg +
                                ", not valid Direction=" + to_string(dir) + "."));
        // Checking the direction has not been defined before.
        AssertThrow (parsed_dirs.find(dir) == parsed_dirs.cend(),
                     ExcMessage("Parsing knot vectors for " + parsing_msg +
                                ", Direction=" + to_string(dir) + " defined"
                                " more than once."));
        parsed_dirs.insert(dir);

        knots[dir] = ke->get_values_vector<Real>();

        const auto size = ke->get_attribute<int>("Size");
        // Checking that the specified size matches with the actual vector size.
        AssertThrow (size == knots[dir].size(),
                     ExcMessage("Parsing knot vectors for " + parsing_msg +
                                ", in Direction=" + to_string(dir) +
                                " Size=" + to_string(size) + " do not match "
                                "with the vector size."));
    }

    const auto grid = GridType::create(knots);
    container->insert_object<GridType>(grid, object_id);
}



template <int dim, int range, int rank>
void
ObjectsContainerParser::
parse_spline_space(const shared_ptr<XMLElement> xml_elem,
                   const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "SplineSpace",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Range") == range,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
    Assert (xml_elem->get_attribute<int>("Rank") == rank,
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

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");

    const auto grid_tag = xml_elem->get_single_element("Grid");
    const auto grid_id = grid_tag->get_attribute<Index>("GetFromIgaObjectId");

    const string parsing_msg = Self_::get_type_id_string("SplineSpace", object_id,
                                                        {{dim, range, rank}});

    // Checking the grid with proper dimension and id exists.
    AssertThrow (container->is_object<GridType> (grid_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching definition" +
                            "for " + Self_::get_type_id_string("Grid", grid_id,
                            SafeSTLVector<Index>(dim)) + "."));
    const auto grid = container->get_object<GridType>(grid_id);
    const auto grid_num_intervals = grid->get_num_intervals();

    // Parsing spline space components.
    const auto comps_elem = xml_elem->get_single_element("SplineSpaceComponents");

    // Checking that the right number of components is defined.
    const auto comp_children = comps_elem->get_children_elements("SplineSpaceComponent");
    AssertThrow (n_components == comp_children.size(),
                 ExcMessage("Parsing " + parsing_msg + ", the number of "
                            "SplineSpaceComponent XML elements does not match "
                            "with the number of components of the space."));

    SafeSTLSet<Index> parsed_comps;
    for (const auto &comp_elem : comp_children)
    {
        const auto comp_id = comp_elem->get_attribute<Index>("ComponentId");
        // Checking the <= 0 comp_id < n_components
        AssertThrow (comp_id >= 0 && comp_id < n_components,
                     ExcMessage("Parsing SplineSpaceComponents for " + parsing_msg +
                                ", not valid ComponentId=" + to_string(comp_id) + "."));
        // Checking the component id has not been defined before.
        AssertThrow (parsed_comps.find(comp_id) == parsed_comps.cend(),
                     ExcMessage("Parsing SplineSpaceComponents for " + parsing_msg +
                                ", ComponentId=" + to_string(comp_id) + " defined"
                                " more than once."));
        parsed_comps.insert(comp_id);

        const auto degree_elem = comp_elem->get_single_element("Degrees");
        const auto degs_vector = degree_elem->get_values_vector<Index>();
        // Check here that degs_vector.size() == dim
        AssertThrow (degs_vector.size() == dim,
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
        AssertThrow (int_mults_elem.size() == dim,
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
            AssertThrow (dir >= 0 && dir < dim,
                         ExcMessage("Parsing InteriorMultiplicities for " +
                                    parsing_msg +
                                    " SplineSpaceComponent ComponentId=" +
                                    to_string(comp_id) +", not valid "
                                    "Direction=" + to_string(dir) + "."));

            // Checking the direction has not been defined before.
            AssertThrow (parsed_dirs.find(dir) == parsed_dirs.cend(),
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
            AssertThrow (size == mults.size(),
                         ExcMessage("Parsing InteriorMultiplicities for " +
                                    parsing_msg +
                                    " SplineSpaceComponent ComponentId=" +
                                    to_string(comp_id) +" in Direction=" +
                                    to_string(dir) + ", Size=" +
                                    to_string(size) + " do not match "
                                    "with the vector size."));

            // Check here that the multiplicities match with the grid.
            AssertThrow ((grid_num_intervals[dir] - 1) == mults.size(),
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
            AssertThrow (periodic_vector.size() == dim,
                 ExcMessage("Parsing " + parsing_msg +
                            ", in SplineSpaceComponent ComponentId=" +
                            to_string(comp_id) + " the number of "
                            " values defined for Periodicity does not "
                            "match with the Space dimension."));

            for (int d = 0; d < dim; ++d)
                period_table[comp_id][d] = periodic_vector[d];
        }
    } // Spline Space components

    const auto spline_space = SplineSpaceType::create
            (deg_table, grid, mult_table, period_table);
    container->insert_object<SplineSpaceType>(spline_space, object_id);
}



template <int dim, int range, int rank>
void
ObjectsContainerParser::
parse_bspline(const shared_ptr<XMLElement> xml_elem,
              const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "BSpline",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Range") == range,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
    Assert (xml_elem->get_attribute<int>("Rank") == rank,
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


    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    const auto ssp_tag = xml_elem->get_single_element("SplineSpace");
    const auto ssp_id = ssp_tag->get_attribute<Index>("GetFromIgaObjectId");

    const string parsing_msg = Self_::get_type_id_string("SplineSpace", object_id,
                                                        {{dim, range, rank}});

    // Checking the spline space with proper dimension and id exists.
    AssertThrow (container->is_object<SplineSpaceType> (ssp_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("SplineSpace", ssp_id,
                            {{dim, range, rank}}) + "."));
    const auto ssp = container->get_object<SplineSpaceType>(ssp_id);
    Assert (ssp != nullptr, ExcNullPtr());

    // Parsing end bevaviour, if exists
    if (xml_elem->has_element("EndBehaviour"))
    {

        // Checking that the right number of components is defined.
        const auto eb_elems = xml_elem->get_single_element("EndBehaviour")
                ->get_children_elements("EndBehaviour");
        AssertThrow (n_components == eb_elems.size(),
                     ExcMessage("Parsing " + parsing_msg + ", the number of "
                                "EndBehaviour elements does not match "
                                "with the number of components of the space."));

        const auto &ssp_periodic_table = ssp->get_periodic_table();

        SafeSTLSet<Index> parsed_comps;
        for (const auto &eb : eb_elems)
        {
            const auto comp_id = eb->get_attribute<Index>("ComponentId");
            // Checking the <= 0 comp_id < n_components
            AssertThrow (comp_id >= 0 && comp_id < n_components,
                         ExcMessage("Parsing EndBehaviour for " + parsing_msg +
                                    ", not valid ComponentId=" + to_string(comp_id) + "."));
            // Checking the component id has not been defined before.
            AssertThrow (parsed_comps.find(comp_id) == parsed_comps.cend(),
                         ExcMessage("Parsing EndBehaviour for " + parsing_msg +
                                    ", ComponentId=" + to_string(comp_id) + " defined"
                                    " more than once."));
            parsed_comps.insert(comp_id);

            const auto string_vec =  eb->get_values_vector<string>();

            // Checking the dimension of the end behaviour vector.
            AssertThrow (string_vec.size() == dim,
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
                    AssertThrow (!ssp_periodic,
                      ExcMessage("Parsing " + parsing_msg + ", in "
                      "EndBehaviour ComponentId=" + to_string(comp_id) +
                      " Direction=" + to_string(dir) + ", behaviour " +
                      sv + " do not match with Spline Space, that is "
                      "periodic for this component and direction."));
                    end_beh_table[comp_id][dir] = BasisEndBehaviour::interpolatory;
                }
                else if (sv == "end_knots")
                {
                    AssertThrow (false,
                      ExcMessage("Parsing " + parsing_msg + ", in "
                      "EndBehaviour ComponentId=" + to_string(comp_id) +
                      " Direction=" + to_string(dir) + ", behaviour " +
                      sv + ". CURRENTLY end_knots cannot be selected from"
                      " the XML input file."));

                    AssertThrow (!ssp_periodic,
                      ExcMessage("Parsing " + parsing_msg + ", in "
                      "EndBehaviour ComponentId=" + to_string(comp_id) +
                      " Direction=" + to_string(dir) + ", behaviour " +
                      sv + " do not match with Spline Space, that is "
                      "periodic for this component and direction."));
                    end_beh_table[comp_id][dir] = BasisEndBehaviour::end_knots;
                }
                else if (sv == "periodic")
                {
                    AssertThrow (ssp_periodic,
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

    const auto bspline = BSplineType::create(ssp, end_beh_table);
    container->insert_object<RefSpaceType>(bspline, object_id);
}



template <int dim>
void
ObjectsContainerParser::
parse_identity_grid_function(const shared_ptr<XMLElement> xml_elem,
                       const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "IdentityGridFunction",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));

    using IdentityGridFunctionType = grid_functions::IdentityGridFunction<dim>;
    using GridFunctionType = GridFunction<dim, dim>;
    using GridType = Grid<dim>;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    const string parsing_msg = Self_::get_type_id_string("IdentityGridFunction",
                               object_id, SafeSTLVector<int>(1, dim));

    const auto gr_tag = xml_elem->get_single_element("Grid");
    const auto gr_id = gr_tag->get_attribute<Index>("GetFromIgaObjectId");

    // Checking the grid with proper dimension and id exists.
    AssertThrow (container->is_object<GridType> (gr_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching definition" +
                            "for " + Self_::get_type_id_string("Grid", gr_id,
                            SafeSTLVector<Index>(dim)) + "."));

    const auto grid = container->get_object<GridType>(gr_id);

    const auto id_func = IdentityGridFunctionType::create(grid);
    container->insert_object<GridFunctionType>(id_func, object_id);
}



template <int dim, int space_dim>
void
ObjectsContainerParser::
parse_constant_grid_function(const shared_ptr<XMLElement> xml_elem,
                       const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "ConstantGridFunction",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Spacedim") == space_dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Spacedim"), space_dim));

    using GridType = Grid<dim>;
    using ConstGridFunctionType = grid_functions::ConstantGridFunction<dim, space_dim>;
    using GridFunctionType = GridFunction<dim, space_dim>;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    const string parsing_msg = Self_::get_type_id_string("IgGridFunction",
                               object_id, {{dim, space_dim}});

    // Gettting grid.
    const auto gr_tag = xml_elem->get_single_element("Grid");
    const auto gr_id = gr_tag->get_attribute<Index>("GetFromIgaObjectId");

    // Checking the grid with proper dimension and id exists.
    AssertThrow (container->is_object<GridType> (gr_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching definition" +
                            "for " + Self_::get_type_id_string("Grid", gr_id,
                            SafeSTLVector<Index>(dim)) + "."));
    const auto grid = container->get_object<GridType>(gr_id);

    // Parsing values.
    const auto vals_tag = xml_elem->get_single_element("Values");

    const auto vals_vec = vals_tag->get_values_vector<Real>();
    AssertThrow (vals_vec.size() == space_dim,
                 ExcMessage("Parsing " + parsing_msg + ", the number of "
                            "components in Values XML does not match "
                            "with the number of components of the GridFunction."));

    typename GridFunctionType::Value values;
    for (int c = 0; c < space_dim; ++c)
        values[c] = vals_vec[c];

    const auto cgf = ConstGridFunctionType::create (grid, values);
    container->insert_object<GridFunctionType>(cgf, object_id);
}



template <int dim, int space_dim>
void
ObjectsContainerParser::
parse_ig_grid_function(const shared_ptr<XMLElement> xml_elem,
                       const bool &first_parsing,
                       const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "IgGridFunction",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Spacedim") == space_dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Spacedim"), space_dim));

    static const int rank = 1;
    using IgGridFunctionType = IgGridFunction<dim, space_dim>;
    using GridFunctionType = GridFunction<dim, space_dim>;
    using RefSpaceType     = ReferenceSpaceBasis<dim, space_dim, rank>;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    const string parsing_msg = Self_::get_type_id_string("IgGridFunction",
                               object_id, {{dim, space_dim}});

    const auto rs_tag = xml_elem->get_single_element("ReferenceSpaceBasis");
    const auto rs_id = rs_tag->get_attribute<Index>("GetFromIgaObjectId");

    // Checking the spline space with proper dimension and id exists.
    if (!container->is_object<RefSpaceType> (rs_id))
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
            AssertThrow (false,
                         ExcMessage("Parsing " + parsing_msg + " not matching "
                                    "definition for " +
                                    Self_::get_type_id_string("ReferenceSpaceBasis", rs_id,
                                    {{dim, space_dim, rank}}) + "."));

    }
    else if (!first_parsing && container->is_object<GridFunctionType> (object_id))
        // We check if the ig grid function was already parsed.
        return;

    const auto rs = container->get_object<RefSpaceType>(rs_id);
    Assert (rs != nullptr, ExcNullPtr());

    const string dofs_property = parse_dofs_property(xml_elem);
    const auto dof_distribution = rs->get_spline_space()->get_dof_distribution();
    AssertThrow (dof_distribution->is_property_defined(dofs_property),
                 ExcMessage("Parsing " + parsing_msg + " dofs property \"" +
                            dofs_property + "\" not defined for " +
                            Self_::get_type_id_string("ReferenceSpaceBasis", rs_id,
                            {{dim, space_dim, rank}}) + "."));

    const auto &global_dofs = dof_distribution->get_global_dofs(dofs_property);
    const auto ig_coefs = parse_ig_coefficients(xml_elem, parsing_msg, global_dofs);

    AssertThrow (rs->get_num_basis() == ig_coefs.size(),
                 ExcMessage("Parsing " + parsing_msg + " the cardinality "
                            "of the ReferenceSpaceBasis (" +
                            to_string(rs->get_num_basis()) + ") is "
                            "different to the dimension of the "
                            "IgCoefficients (" + to_string(ig_coefs.size())
                            + ")."));

    const auto igf = IgGridFunctionType::create(rs, ig_coefs, dofs_property);
    container->insert_object<GridFunctionType>(igf, object_id);
}



template <int dim, int range, int rank>
void
ObjectsContainerParser::
parse_nurbs(const shared_ptr<XMLElement> xml_elem,
            const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "NURBS",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Range") == range,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
    Assert (xml_elem->get_attribute<int>("Rank") == rank,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));

    using BSplineType = BSpline<dim, range, rank>;
    using NURBSType   = NURBS<dim, range, rank>;
    using WeightIgFunctionType = IgGridFunction<dim, 1>;
    using WeightFunctionType = GridFunction<dim, 1>;
    using RefSpaceType = ReferenceSpaceBasis<dim, range, rank>;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");

    const string parsing_msg = Self_::get_type_id_string("NURBS",
                               object_id, {{dim, range, rank}});

    // Parsing BSpline
    const auto bs_tag = xml_elem->get_single_element("BSpline");
    const auto bs_id = bs_tag->get_attribute<Index>("GetFromIgaObjectId");
    // Checking the bspline with the proper dimensions and id exists.
    AssertThrow (container->is_object<RefSpaceType> (bs_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("BSpline", bs_id,
                            {{dim, range, rank}}) + "."));

    const auto rs = container->get_object<RefSpaceType>(bs_id);
    Assert (rs != nullptr, ExcNullPtr());
    // Checking that the reference space is a BSpline.
    const auto bs = dynamic_pointer_cast<BSplineType>(rs);
    Assert (bs != nullptr, ExcNullPtr());
    AssertThrow (rs->is_bspline(),
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("BSplineSpaceBasis", bs_id,
                            {{dim, range, rank}}) + ". It is a "
                            "ReferenceBasisSpace, but not a BSpline."
                            "BSpline space basis must be used for "
                            "defining IgGridFunctions used as weight "
                            "functions for NURBS space basis."));

    // Parsing Weight function
    const auto wg_tag = xml_elem->get_single_element("WeightFunction");
    const auto wg_id = wg_tag->get_attribute<Index>("GetFromIgaObjectId");
    // Checking that the grid function with the proper dimensions exists.
    AssertThrow (container->is_object<WeightFunctionType> (wg_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("IgGridFunction", wg_id,
                            {{dim, 1}}) +  ". The error could be caused "
                            "by the fact the object does not correspond "
                            "to IgGridFunction built upon a BSpline. "
                            "Currently igatools only allows to build "
                            "weight functions based on scalar BSpline "
                            "space basis."));

    const auto gf = container->get_object<WeightFunctionType>(wg_id);
    Assert (gf != nullptr, ExcNullPtr());
    const auto wf = dynamic_pointer_cast<WeightIgFunctionType>(gf);
    // Checking that the grid function is a ig grid function.
    AssertThrow (wf != nullptr,
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("IgGridFunction", wg_id,
                            {{dim, 1}}) + ". It is a"
                            "GridFunction, but not a IgGridFunction"));

    // Checking that the weight function was defined upon a BSpline.
    const auto w_func_basis = wf->get_basis();
    Assert(wf->get_basis()->is_bspline(),
           ExcMessage("Parsing " + parsing_msg + ", " +
                      Self_::get_type_id_string("IgGridFunction", wg_id,
                      {{dim, 1}}) + ". It is based on a NURBS space basis"
                      ", but must be based on BSpline space basis."));

    // Checking that the grids of the bspline and weight function match.
    AssertThrow(bs->get_grid() == wf->get_grid(),
                ExcMessage("Parsing " + parsing_msg + ", mismatching "
                           "grids for the BSpline and WeightFunction."));

    // Checking that the number of basis functions and weights match.
    const auto &n_basis_table = bs->get_spline_space()->get_dof_distribution()->get_num_dofs_table();
    int comp_id = 0;
    for (const auto &n_basis_comp : n_basis_table)
    {
        AssertThrow(n_basis_comp ==
                    w_func_basis->get_spline_space()->get_dof_distribution()->get_num_dofs_table()[0],
               ExcMessage("Parsing " + parsing_msg + ", mismatching number of basis functions and weight "
                          "coefficients for scalar component " + to_string(comp_id++)));
    }

    const auto nurbs = NURBSType::create(bs, wf);
    Assert (nurbs != nullptr, ExcNullPtr());
    container->insert_object<RefSpaceType>(nurbs, object_id);
}



template <int dim, int codim>
void
ObjectsContainerParser::
parse_domain(const shared_ptr<XMLElement> xml_elem,
             const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "Domain",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Codim") == codim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Codim"), codim));

    static const int space_dim = dim + codim;
    using DomainType = Domain<dim, codim>;
    using GridFuncType = GridFunction<dim, space_dim>;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");

    const string parsing_msg = Self_::get_type_id_string("Domain",
                               object_id, {{dim, codim}});

    const auto gf_tag = xml_elem->get_single_element("GridFunction");
    const auto gf_id = gf_tag->get_attribute<Index>("GetFromIgaObjectId");
    // Checking if the grid function exists.
    AssertThrow (container->is_object<GridFuncType> (gf_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("GridFunction", gf_id,
                            {{dim, space_dim}}) + "."));

    const auto gf = container->get_object<GridFuncType>(gf_id);
    Assert (gf != nullptr, ExcNullPtr());

    const auto name = parse_name (xml_elem);
    const auto domain = DomainType::create(gf, name);

    container->insert_object<DomainType>(domain, object_id);
}



template <int dim, int codim, int range, int rank>
void
ObjectsContainerParser::
parse_ig_function(const shared_ptr<XMLElement> xml_elem,
                  const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "IgFunction",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Range") == range,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
    Assert (xml_elem->get_attribute<int>("Rank") == rank,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));
    Assert (xml_elem->get_attribute<int>("Codim") == codim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Codim"), codim));

    using PhysSpaceType = PhysicalSpaceBasis<dim, range, rank, codim>;
    using IgFunctionType = IgFunction<dim, codim, range, rank>;
    using FunctionType = Function<dim, codim, range, rank>;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");

    const string parsing_msg = Self_::get_type_id_string("IgFunction",
                               object_id, {{dim, codim, range, rank}});

    const auto ps_tag = xml_elem->get_single_element("PhysicalSpaceBasis");
    const auto ps_id = ps_tag->get_attribute<Index>("GetFromIgaObjectId");

    // Checking if the physical space exists.
    AssertThrow (container->is_object<PhysSpaceType> (ps_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("PhysicalSpaceBasis", ps_id,
                            {{dim, range, rank, codim}}) + "."));
    const auto ps = container->get_object<PhysSpaceType>(ps_id);
    Assert (ps != nullptr, ExcNullPtr());

    const auto name = parse_name (xml_elem);

    const string dofs_property = parse_dofs_property(xml_elem);
    const auto dof_distribution = ps->get_spline_space()->get_dof_distribution();
    AssertThrow (dof_distribution->is_property_defined(dofs_property),
                 ExcMessage("Parsing " + parsing_msg + " dofs property \"" +
                            dofs_property + "\" not defined for " +
                            Self_::get_type_id_string("PhysicalSpaceBasis", ps_id,
                            {{dim, range, rank, codim}}) + "."));

    const auto &global_dofs = dof_distribution->get_global_dofs(dofs_property);
    const auto ig_coefs = parse_ig_coefficients(xml_elem, parsing_msg, global_dofs);

    AssertThrow (ps->get_num_basis() == ig_coefs.size(),
                 ExcMessage("Parsing " + parsing_msg + " the cardinality "
                            "of the PhysicalSpaceBasis (" +
                            to_string(ps->get_num_basis()) + ") is "
                            "different to the dimension of the "
                            "IgCoefficients (" + to_string(ig_coefs.size())
                            + ")."));

    const auto igf = IgFunctionType::create(ps, ig_coefs, dofs_property, name);
    container->insert_object<FunctionType>(igf, object_id);
}



template <int dim, int codim, int range, int rank>
void
ObjectsContainerParser::
parse_constant_function(const shared_ptr<XMLElement> xml_elem,
                        const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "ConstantFunction",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Range") == range,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
    Assert (xml_elem->get_attribute<int>("Rank") == rank,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));
    Assert (xml_elem->get_attribute<int>("Codim") == codim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Codim"), codim));

    using DomainType = Domain<dim, codim>;
    using FunctionType = Function<dim, codim, range, rank>;
    using ConstFunctionType = functions::ConstantFunction<dim, codim, range, rank>;
    using Values = typename FunctionType::Value;
    static const int n_components = Values::size;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");

    const string parsing_msg = Self_::get_type_id_string("ConstantFunction",
                               object_id, {{dim, codim, range, rank}});

    const auto dm_tag = xml_elem->get_single_element("Domain");
    const auto dm_id = dm_tag->get_attribute<Index>("GetFromIgaObjectId");

    AssertThrow (container->is_object<DomainType> (dm_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("Domain", dm_id,
                            {{dim, codim}}) + "."));

    const auto dm = container->get_object<DomainType>(dm_id);
    Assert (dm != nullptr, ExcNullPtr());

    // Parsing values.
    const auto vals_tag = xml_elem->get_single_element("Values");

    const auto vals_vec = vals_tag->get_values_vector<Real>();
    AssertThrow (vals_vec.size() == n_components,
                 ExcMessage("Parsing " + parsing_msg + ", the number of "
                            "components in Values XML does not match "
                            "with the number of components of the Function."));

    typename FunctionType::Value values;
    for (int c = 0; c < n_components; ++c)
        values[c] = vals_vec[c];

    const auto name = parse_name (xml_elem);

    const auto cgf = ConstFunctionType::create (dm, values, name);
    container->insert_object<FunctionType>(cgf, object_id);
}



template <int dim, int codim, int range, int rank>
void
ObjectsContainerParser::
parse_phys_space(const shared_ptr<XMLElement> xml_elem,
                 const shared_ptr<ObjectsContainer> container)
{
    Assert (xml_elem->get_name() == "PhysicalSpaceBasis",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Range") == range,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Range"), range));
    Assert (xml_elem->get_attribute<int>("Rank") == rank,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Rank"), rank));
    Assert (xml_elem->get_attribute<int>("Codim") == codim,
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

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");

    const auto parsing_msg = Self_::get_type_id_string("PhysicalSpaceBasis",
            object_id, {{dim, range, rank}});

    const auto rs_tag = xml_elem->get_single_element("ReferenceSpaceBasis");
    const auto rs_id = rs_tag->get_attribute<Index>("GetFromIgaObjectId");

    AssertThrow (container->is_object<RefSpaceType> (rs_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("ReferenceSpaceBasis", rs_id,
                            {{dim, range, rank}}) + "."));

    const auto rs = container->get_object<RefSpaceType>(rs_id);
    Assert (rs != nullptr, ExcNullPtr());

    const auto dm_tag = xml_elem->get_single_element("Domain");
    const auto dm_id = dm_tag->get_attribute<Index>("GetFromIgaObjectId");

    AssertThrow (container->is_object<DomainType> (dm_id),
                 ExcMessage("Parsing " + parsing_msg + " not matching "
                            "definition for " +
                            Self_::get_type_id_string("Domain", dm_id,
                            {{dim, codim}}) + "."));

    const auto dm = container->get_object<DomainType>(dm_id);
    Assert (dm != nullptr, ExcNullPtr());

    // Checking that the reference space and domain grids match.
    AssertThrow(rs->get_grid() == dm->get_grid_function()->get_grid(),
                ExcMessage("Parsing " + parsing_msg + ", mismatching "
                           "grids for the ReferenceSpaceBasis and Domain."));

    const auto ps = PhysSpaceType::create(rs, dm, transf);
    container->insert_object<PhysSpaceType>(ps, object_id);
}



string
ObjectsContainerParser::
parse_name(const shared_ptr<XMLElement> xml_elem)
{
    if (xml_elem->has_element("Name"))
        return xml_elem->get_single_element("Name")->get_value<string>();
    else
        return string("");
}



string
ObjectsContainerParser::
parse_dofs_property(const shared_ptr<XMLElement> xml_elem)
{
    if (xml_elem->has_element("DofsProperty"))
        return xml_elem->get_single_element("DofsProperty")->get_value<string>();
    else
        return DofProperties::active;
}



string
ObjectsContainerParser::
get_type_id_string(const string &object_type,
                   const Index &object_id,
                   const SafeSTLVector<int> &dims)
{

    string dims_str = object_type + "<";
    for (const auto &d : dims)
        dims_str += to_string(d) + ", ";
    return dims_str.substr(0, dims_str.size() - 2) + ">" +
            " (IgaObjectId=" + to_string(object_id) + ")";
}



IgCoefficients
ObjectsContainerParser::
parse_ig_coefficients(const shared_ptr<XMLElement> xml_elem,
                      const string &parsing_msg,
                      const set<Index> &space_global_dofs)
{
    Assert (xml_elem->has_element("IgCoefficients"),
            ExcMessage("IgCoefficients XML element not present."));

    const auto ig_elem = xml_elem->get_single_element("IgCoefficients");
    const auto size = ig_elem->get_attribute<Index>("Size");

    // Gettint the indices.
    const auto ig_ind_elem = ig_elem->get_single_element("Indices");
    const auto ig_coefs_ind_vec = ig_ind_elem->get_values_vector<Index>();
    AssertThrow (size == ig_coefs_ind_vec.size(),
                 ExcMessage("Parsing IgCoefficients indices for " + parsing_msg +
                            ", Size=" + to_string(size) + " do not match "
                            "with the vector size."));

    // Gettint the values.
    const auto ig_val_elem = ig_elem->get_single_element("Values");
    const auto ig_coefs_val_vec = ig_val_elem->get_values_vector<Real>();
    AssertThrow (size == ig_coefs_val_vec.size(),
                 ExcMessage("Parsing IgCoefficients values for " + parsing_msg +
                            ", Size=" + to_string(size) + " do not match "
                            "with the vector size."));

    set<Index> indices;
    for (const auto &i : ig_coefs_ind_vec)
        indices.insert(i);
    copy(ig_coefs_ind_vec.cbegin(), ig_coefs_ind_vec.cend(),
              inserter(indices, indices.begin()));

    // Checking that there are not repeated indices.
    AssertThrow (size == indices.size(),
                 ExcMessage("Parsing IgCoefficients for " + parsing_msg +
                            ", not valid indices vector parsed. Repeated "
                            "indices may found."));

    IgCoefficients ig_coefs (indices);

    // Checking if the parsed indices match with space_global_dofs and
    // filling the values of the ig coefficients vector.
    const auto end_dofs = space_global_dofs.cend();
    auto ind_it = ig_coefs_ind_vec.cbegin();
    auto val_it = ig_coefs_val_vec.cbegin();
    auto val_end = ig_coefs_val_vec.cend();
    for (; val_it != val_end; ++val_it, ++ind_it)
    {
        ig_coefs[*ind_it] = *val_it;
        AssertThrow (space_global_dofs.find(*ind_it) != end_dofs,
                 ExcMessage("Parsing IgCoefficients for " + parsing_msg +
                            ", " + to_string(*ind_it) + " is not a valid index."));
    }

    return ig_coefs;
}

IGA_NAMESPACE_CLOSE

#endif // XML_IO
