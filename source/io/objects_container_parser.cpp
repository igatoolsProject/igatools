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

#include <igatools/geometry/grid.h>
#include <igatools/basis_functions/spline_space.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/functions/function.h>
#include <igatools/functions/ig_function.h>

#include <boost/fusion/algorithm/iteration/for_each.hpp>

using std::string;
using std::to_string;
using std::shared_ptr;

IGA_NAMESPACE_OPEN


ObjectsContainerParser::
ObjectsContainerParser(const string &file_path)
  :
  file_parser_(XMLFileParser::create(file_path))
{}



auto
ObjectsContainerParser::
create(const string &file_path) -> SelfPtr_
{
    return SelfPtr_ (new Self_ (file_path));
}



shared_ptr<ObjectsContainer>
ObjectsContainerParser::
parse(const string &schema_file) const
{
    const shared_ptr<XMLElement> xml_elem = file_parser_->parse(schema_file);
    const auto container = ObjectsContainer::create();

    // Checking for repeated iga object ids.
    SafeSTLSet<Index> object_ids;
    for (const auto &el : xml_elem->get_children_elements())
    {
        const Index obj_id = el->get_attribute<Index>("IgaObjectId");
        AssertThrow (object_ids.find(obj_id) == object_ids.cend(),
                ExcMessage("IgaObjectId " + to_string(obj_id) + " is "
                           "used more than once."));
        object_ids.insert(obj_id);
    }


    this->parse_grids (xml_elem, container);
    this->parse_spline_spaces (xml_elem, container);
    this->parse_bsplines (xml_elem, container);
    this->parse_grid_functions_and_nurbs (xml_elem, container);
    this->parse_domains (xml_elem, container);
    this->parse_phys_spaces (xml_elem, container);
    this->parse_ig_functions (xml_elem, container);
    return container;
}



void
ObjectsContainerParser::
parse_grids(const shared_ptr<XMLElement> xml_elem,
            const shared_ptr<ObjectsContainer> container) const
{
    const auto grid_elems = xml_elem->get_children_elements("Grid");
    for (const auto &ge : grid_elems)
    {
        const int grid_dim = ge->get_attribute<int>("Dim");

        using ValidGridPtrs = typename ObjectsContainer::ValidGridPtrs;
        ValidGridPtrs valid_grid_ptr_types;

        bool found = false;
        boost::fusion::for_each(valid_grid_ptr_types, [&](const auto &grid_ptr_type)
        {
            if (found)
                return;

            using GridType = typename std::remove_reference<decltype(grid_ptr_type)>::type::element_type;
            static const int dim = GridType::dim;

            if (grid_dim == dim)
            {
                found = true;
                parse_grid<dim>(ge, container);
            }
        });

        // Grid dimension not found
        AssertThrow (found,
                     ExcMessage("Parsing Grid with IgaObjectId=" +
                                to_string(ge->get_attribute<Index>("IgaObjectId")) +
                                " not valid dimension "
                                + to_string(ge->get_attribute<int>("Dim")) +
                                " for a grid. Maybe igatools was not "
                                "instantiated for this dimension."));
    }
}



template <int dim>
void
ObjectsContainerParser::
parse_grid(const shared_ptr<XMLElement> xml_elem,
           const std::shared_ptr<ObjectsContainer> container) const
{
    Assert (xml_elem->get_name() == "Grid",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    AssertThrow (!container->is_id_present(object_id),
                 ExcMessage("Parsing Grid already defined IgaObjectId=" +
                            to_string(object_id) + "."));

    // Check here that the object id is not already defined.

    const auto knots_children = xml_elem->get_children_elements("Knots");
    AssertThrow (dim == knots_children.size(),
                 ExcMessage("Parsing Grid with IgaObjectId=" +
                            to_string(object_id) + " the number of knots vectors does " +
                            "match with Dim of the Grid."));

    using GridType = Grid<dim>;

    SafeSTLArray<SafeSTLVector<Real>, dim> knots;

    SafeSTLSet<Index> parsed_dirs;

    for (const auto &ke : knots_children)
    {
        const auto dir = ke->get_attribute<Index>("Direction");
        AssertThrow (parsed_dirs.find(dir) == parsed_dirs.cend(),
                     ExcMessage("Parsing Grid with IgaObjectId=" +
                                to_string(object_id) + " knot vector for "
                                "Direction " + to_string(dir) + " defined"
                                " more than once."));
        // TODO: to check here that the direction is lower than dir.
        parsed_dirs.insert(dir);

        knots[dir] = ke->get_values_vector<Real>();

        AssertThrow (ke->get_attribute<int>("Size") == knots[dir].size(),
                     ExcMessage("Parsing Grid with IgaObjectId=" +
                                to_string(object_id) + " knot vector Size for Direction " +
                                to_string(dir) + " does not match with the specified one."));
    }

    const auto grid = GridType::create(knots);

    container->insert_object<GridType>(grid, object_id);
}



void
ObjectsContainerParser::
parse_spline_spaces(const shared_ptr<XMLElement> xml_elem,
                    const shared_ptr<ObjectsContainer> container) const
{
    for (const auto &ssp : xml_elem->get_children_elements("SplineSpace"))
    {
        const int ssp_dim = ssp->get_attribute<int>("Dim");
        const int ssp_range = ssp->get_attribute<int>("Range");
        const int ssp_rank = ssp->get_attribute<int>("Rank");

        using ValidSplineSpacePtrs = typename ObjectsContainer::ValidSplineSpacePtrs;
        ValidSplineSpacePtrs valid_ssp_ptr_types;

        bool found = false;
        boost::fusion::for_each(valid_ssp_ptr_types, [&](const auto &ssp_ptr_type)
        {
            if (found)
                return;

            using SplineSpaceType = typename
                    std::remove_reference<decltype(ssp_ptr_type)>::type::element_type;
            static const int dim = SplineSpaceType::dim;
            static const int range = SplineSpaceType::range;
            static const int rank = SplineSpaceType::rank;

            if (ssp_dim == dim && ssp_range == range && ssp_rank == rank)
            {
                found = true;
                parse_spline_space<dim, range, rank>(ssp, container);
            }
        });
    }
}



template <int dim, int range, int rank>
void
ObjectsContainerParser::
parse_spline_space(const shared_ptr<XMLElement> xml_elem,
                   const std::shared_ptr<ObjectsContainer> container) const
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

    const bool default_periodicity = false;

    DegreeTable deg_table;
    MultiplicityTable mult_table;
    PeriodicityTable period_table;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    AssertThrow (!container->is_id_present(object_id),
                 ExcMessage("Parsing SplineSpace already defined IgaObjectId=" +
                            to_string(object_id) + "."));

    const auto grid_tag = xml_elem->get_single_element("Grid");
    const auto grid_id = grid_tag->get_attribute<Index>("GetFromIgaObjectId");

    AssertThrow (container->is_object<GridType> (grid_id),
                 ExcMessage("Parsing SplineSpace with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a grid with the expected dimension."));

    const auto grid = container->get_object<GridType>(grid_id);

    const auto comps_elem = xml_elem->get_single_element("SplineSpaceComponents");
    // Check that comps_elem == n_components
    for (const auto &comp_elem : comps_elem->get_children_elements("SplineSpaceComponent"))
    {
        const auto comp_id = comp_elem->get_attribute<Index>("ComponentId");
        // To check here 0 <= comp_id < n_components and is not repeated.

        const auto degree_elem = comp_elem->get_single_element("Degrees");
        const auto degs_vector = degree_elem->get_values_vector<Index>();
        // Check here that degs_vector.size() == dim
        for (int d = 0; d < dim; ++d)
            deg_table[comp_id][d] = degs_vector[d];

        const auto int_mults_elem =
                comp_elem->get_single_element("InteriorMultiplicities")
                ->get_children_elements("InteriorMultiplicities");
        // Check here that  int_muls_elem.size() == dim

        for (const auto &im : int_mults_elem)
        {
            const auto dir = im->get_attribute<Index>("Direction");
//            AssertThrow (parsed_dirs.find(dir) == parsed_dirs.cend(),
//                         ExcMessage("Parsing Grid with IgaObjectId=" +
//                                    to_string(object_id) + " knot vector for "
//                                    "Direction " + to_string(dir) + " defined"
//                                    " more than once."));
//            parsed_dirs.insert(dir);

            const auto mults = im->get_values_vector<Index>();
            mult_table[comp_id].copy_data_direction(dir, mults);
            // Check here that the multiplicities match with the grid.

//            AssertThrow (im->get_attribute<int>("Size") == mults.size(),
//                         ExcMessage("Parsing Grid with IgaObjectId=" +
//                                    to_string(object_id) + " knot vector Size for Direction " +
//                                    to_string(dir) + " does not match with the specified one."));
        }

        if (comp_elem->has_element("Periodicity"))
        {
            const auto periodic_vector = comp_elem->
                    get_single_element("Periodicity")->get_values_vector<bool>();
            // Check dimension
            for (int d = 0; d < dim; ++d)
                period_table[comp_id][d] = periodic_vector[d];
        }
        else
        {
            for (int d = 0; d < dim; ++d)
                period_table[comp_id][d] = default_periodicity;

        }
    } // Spline Space components


    const auto spline_space = SplineSpaceType::create (deg_table, grid, mult_table, period_table);

    container->insert_object<SplineSpaceType>(spline_space, object_id);
}



void
ObjectsContainerParser::
parse_bsplines(const shared_ptr<XMLElement> xml_elem,
               const shared_ptr<ObjectsContainer> container) const
{
    for (const auto &bs : xml_elem->get_children_elements("BSpline"))
    {
        const int bs_dim = bs->get_attribute<int>("Dim");
        const int bs_range = bs->get_attribute<int>("Range");
        const int bs_rank = bs->get_attribute<int>("Rank");

        using ValidBSplinePtrs = typename ObjectsContainer::ValidSplineSpacePtrs;
        ValidBSplinePtrs valid_bs_ptr_types;

        bool found = false;
        boost::fusion::for_each(valid_bs_ptr_types, [&](const auto &bs_ptr_type)
        {
            if (found)
                return;

            using BSplineType = typename
                    std::remove_reference<decltype(bs_ptr_type)>::type::element_type;
            static const int dim = BSplineType::dim;
            static const int range = BSplineType::range;
            static const int rank = BSplineType::rank;

            if (bs_dim == dim && bs_range == range && bs_rank == rank)
            {
                found = true;
                parse_bspline<dim, range, rank>(bs, container);
            }
        });
    }
}



template <int dim, int range, int rank>
void
ObjectsContainerParser::
parse_bspline(const shared_ptr<XMLElement> xml_elem,
              const std::shared_ptr<ObjectsContainer> container) const
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
    using EndBehaviour = typename SplineSpaceType::EndBehaviour;
    static const int n_components = SplineSpaceType::n_components;
    EndBehaviourTable end_beh_table;

    // Initializing to default values.
    for (auto &eb_c : end_beh_table)
        for (auto &eb : eb_c)
            eb = BasisEndBehaviour::interpolatory;


    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    AssertThrow (!container->is_id_present(object_id),
                 ExcMessage("Parsing BSplineSpace already defined IgaObjectId=" +
                            to_string(object_id) + "."));

    const auto ssp_tag = xml_elem->get_single_element("SplineSpace");
    const auto ssp_id = ssp_tag->get_attribute<Index>("GetFromIgaObjectId");
    AssertThrow (container->is_object<SplineSpaceType> (ssp_id),
                 ExcMessage("Parsing BSpline with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a SplineSpace with the expected dimension."));

    const auto ssp = container->get_object<SplineSpaceType>(ssp_id);
    Assert (ssp != nullptr, ExcNullPtr());

    if (xml_elem->has_element("EndBehaviour"))
    {
        const auto eb_elems = xml_elem->get_single_element("EndBehaviour")
                ->get_children_elements("EndBehaviour");
        // Check eb_elems.size() == n_components

        for (const auto &eb : eb_elems)
        {
            const auto comp_id = eb->get_attribute<Index>("ComponentId");
            // Check component
            const auto string_vec =  eb->get_values_vector<string>();
            // check string_vec.size()
            Index d = 0;
            for (const auto &sv : string_vec)
            {
                if (sv == "interpolatory")
                {
                    // Check spline space periodicity
                    end_beh_table[comp_id][d] = BasisEndBehaviour::interpolatory;
                }
                else if (sv == "end_knots")
                {
                    // Check spline space periodicity
                    end_beh_table[comp_id][d] = BasisEndBehaviour::end_knots;
                }
                else if (sv == "periodic")
                {
                    // Check spline space periodicity
                    end_beh_table[comp_id][d] = BasisEndBehaviour::periodic;
                }
                else
                {
                    // throw error
                }
                ++d;
            }
        }
    }
    const auto bspline = BSplineType::create(ssp, end_beh_table);

    container->insert_object<RefSpaceType>(bspline, object_id);
}



void
ObjectsContainerParser::
parse_grid_functions_and_nurbs(const shared_ptr<XMLElement> xml_elem,
                               const shared_ptr<ObjectsContainer> container) const
{
    AssertThrow (false, ExcNotImplemented());
    // Call properly here parse ig_grid_function and parse_nurbs
}

template <int dim, int range>
void
ObjectsContainerParser::
parse_ig_grid_function(const shared_ptr<XMLElement> xml_elem,
                       const std::shared_ptr<ObjectsContainer> container) const
{
    AssertThrow (false, ExcNotImplemented());
}



template <int dim, int range, int rank>
void
ObjectsContainerParser::
parse_nurbs(const shared_ptr<XMLElement> xml_elem,
            const std::shared_ptr<ObjectsContainer> container) const
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
    AssertThrow (!container->is_id_present(object_id),
                 ExcMessage("Parsing NURBS already defined IgaObjectId=" +
                            to_string(object_id) + "."));

    // Parsing BSpline
    const auto bs_tag = xml_elem->get_single_element("BSpline");
    const auto bs_id = bs_tag->get_attribute<Index>("GetFromIgaObjectId");
    AssertThrow (container->is_object<RefSpaceType> (bs_id),
                 ExcMessage("Parsing NURBS with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a BSpline with the expected dimension."));

    const auto rs = container->get_object<RefSpaceType>(bs_id);
    Assert (rs != nullptr, ExcNullPtr());
    const auto bs = std::dynamic_pointer_cast<BSplineType>(rs);
    AssertThrow (bs != nullptr,
                 ExcMessage("Parsing NURBS with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a BSpline with the expected dimension."));

    // Parsing Weight function
    const auto wg_tag = xml_elem->get_single_element("WeightFunction");
    const auto wg_id = bs_tag->get_attribute<Index>("GetFromIgaObjectId");
    AssertThrow (container->is_object<WeightFunctionType> (wg_id),
                 ExcMessage("Parsing NURBS with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a GridFunction with the expected dimension."));

    const auto gf = container->get_object<WeightFunctionType>(wg_id);
    Assert (gf != nullptr, ExcNullPtr());
    const auto wf = std::dynamic_pointer_cast<WeightIgFunctionType>(gf);
    AssertThrow (wf != nullptr,
                 ExcMessage("Parsing NURBS with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a IgGridFunction with the expected dimension."));

    const auto nurbs = NURBSType::create(bs, wf);
    Assert (nurbs != nullptr, ExcNullPtr());
    container->insert_object<RefSpaceType>(nurbs, object_id);
}



void
ObjectsContainerParser::
parse_domains(const shared_ptr<XMLElement> xml_elem,
              const shared_ptr<ObjectsContainer> container) const
{
    for (const auto &dm : xml_elem->get_children_elements("Domain"))
    {
        const int dm_dim = dm->get_attribute<int>("Dim");
        const int dm_codim = dm->get_attribute<int>("Codim");

        using ValidDomainPtrs = typename ObjectsContainer::ValidDomainPtrs;
        ValidDomainPtrs valid_dm_ptr_types;

        bool found = false;
        boost::fusion::for_each(valid_dm_ptr_types, [&](const auto &dm_ptr_type)
        {
            if (found)
                return;

            using DomainType = typename
                    std::remove_reference<decltype(dm_ptr_type)>::type::element_type;
            static const int dim = DomainType::dim;
            static const int space_dim = DomainType::space_dim;
            static const int codim = space_dim - dim;

            if (dm_dim == dim && dm_codim == codim)
            {
                found = true;
                parse_domain<dim, codim>(dm, container);
            }
        });
    }
}



void
ObjectsContainerParser::
parse_phys_spaces(const shared_ptr<XMLElement> xml_elem,
                  const shared_ptr<ObjectsContainer> container) const
{
    for (const auto &ps : xml_elem->get_children_elements("PhysicalSpaceBasis"))
    {
        const int ps_dim = ps->get_attribute<int>("Dim");
        const int ps_codim = ps->get_attribute<int>("Codim");
        const int ps_range = ps->get_attribute<int>("Range");
        const int ps_rank = ps->get_attribute<int>("Rank");

        using ValidPhysSpacePtrs = typename ObjectsContainer::ValidPhysSpacePtrs;
        ValidPhysSpacePtrs valid_ps_ptr_types;

        bool found = false;
        boost::fusion::for_each(valid_ps_ptr_types, [&](const auto &ps_ptr_type)
        {
            if (found)
                return;

            using PhysSpaceType = typename
                    std::remove_reference<decltype(ps_ptr_type)>::type::element_type;
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
    }
}



void
ObjectsContainerParser::
parse_ig_functions(const shared_ptr<XMLElement> xml_elem,
                  const shared_ptr<ObjectsContainer> container) const
{
    for (const auto &fn : xml_elem->get_children_elements("IgFunction"))
    {
        const int fn_dim = fn->get_attribute<int>("Dim");
        const int fn_codim = fn->get_attribute<int>("Codim");
        const int fn_range = fn->get_attribute<int>("Range");
        const int fn_rank = fn->get_attribute<int>("Rank");

        using ValidFunctionPtrs = typename ObjectsContainer::ValidFunctionPtrs;
        ValidFunctionPtrs valid_fn_ptr_types;

        bool found = false;
        boost::fusion::for_each(valid_fn_ptr_types, [&](const auto &fn_ptr_type)
        {
            if (found)
                return;

            using FunctionType = typename
                    std::remove_reference<decltype(fn_ptr_type)>::type::element_type;
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
    }
}



template <int dim, int codim>
void
ObjectsContainerParser::
parse_domain(const shared_ptr<XMLElement> xml_elem,
             const std::shared_ptr<ObjectsContainer> container) const
{
    Assert (xml_elem->get_name() == "Domain",
            ExcMessage("Invalid XML tag."));

    Assert (xml_elem->get_attribute<int>("Dim") == dim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Dim"), dim));
    Assert (xml_elem->get_attribute<int>("Codim") == codim,
            ExcDimensionMismatch(xml_elem->get_attribute<int>("Codim"), codim));

    using DomainType = Domain<dim, codim>;
    static const int space_dim = dim + codim;
    using GridFuncType = GridFunction<dim, space_dim>;

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    AssertThrow (!container->is_id_present(object_id),
                 ExcMessage("Parsing Domain already defined IgaObjectId=" +
                            to_string(object_id) + "."));

    const auto gf_tag = xml_elem->get_single_element("GridFunction");
    const auto gf_id = gf_tag->get_attribute<Index>("GetFromIgaObjectId");
    AssertThrow (container->is_object<GridFuncType> (gf_id),
                 ExcMessage("Parsing Domain with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a GridFunc with the expected dimension."));

    const auto gf = container->get_object<GridFuncType>(gf_id);
    Assert (gf != nullptr, ExcNullPtr());

    string name = "";
    if (xml_elem->has_element("Name"))
    {
        const auto nm_elem = xml_elem->get_single_element("Name");
        name = nm_elem->get_value<string>();
    }
    const auto domain = DomainType::create(gf, name);

    container->insert_object<DomainType>(domain, object_id);
}



template <int dim, int codim, int range, int rank>
void
ObjectsContainerParser::
parse_ig_function(const shared_ptr<XMLElement> xml_elem,
                  const std::shared_ptr<ObjectsContainer> container) const
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
    AssertThrow (!container->is_id_present(object_id),
                 ExcMessage("Parsing IgFunction already defined IgaObjectId=" +
                            to_string(object_id) + "."));

    const auto ps_tag = xml_elem->get_single_element("PhysicalSpaceBasis");
    const auto ps_id = ps_tag->get_attribute<Index>("GetFromIgaObjectId");
    AssertThrow (container->is_object<PhysSpaceType> (ps_id),
                 ExcMessage("Parsing Ig Function with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a Physical Space Basis with the expected dimension."));
    const auto ps = container->get_object<PhysSpaceType>(ps_id);
    Assert (ps != nullptr, ExcNullPtr());

    string dofs_property = "active";
    if (xml_elem->has_element("DofsProperty"))
    {
        const auto dp_elem = xml_elem->get_single_element("DofsProperty");
        dofs_property = dp_elem->get_value<string>();
    }

    const auto ig_elem = xml_elem->get_single_element("IgCoefficients");
    const auto size = ig_elem->get_attribute<Index>("Size");
    const auto ig_coefs_vec = ig_elem->get_values_vector<Real>();
    // Check ig_coefs_vec.size() == size

    std::set<Index> indices;
    for (int i = 0; i < size; ++i)
        indices.insert(i);
    IgCoefficients ig_coefs (indices);
    for (int i = 0; i < size; ++i)
        ig_coefs[i] = ig_coefs_vec[i];

    string name = "";
    if (xml_elem->has_element("Name"))
    {
        const auto nm_elem = xml_elem->get_single_element("Name");
        name = nm_elem->get_value<string>();
    }

    const auto igf = IgFunctionType::create(ps, ig_coefs, dofs_property, name);
    container->insert_object<FunctionType>(igf, object_id);
}



template <int dim, int codim, int range, int rank>
void
ObjectsContainerParser::
parse_phys_space(const shared_ptr<XMLElement> xml_elem,
                 const std::shared_ptr<ObjectsContainer> container) const
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
        else
        {
            // Throw error
        }
    }

    const auto object_id = xml_elem->get_attribute<Index>("IgaObjectId");
    AssertThrow (!container->is_id_present(object_id),
                 ExcMessage("Parsing Physical Space Basis already defined IgaObjectId=" +
                            to_string(object_id) + "."));

    const auto rs_tag = xml_elem->get_single_element("ReferenceSpaceBasis");
    const auto rs_id = rs_tag->get_attribute<Index>("GetFromIgaObjectId");
    AssertThrow (container->is_object<RefSpaceType> (rs_id),
                 ExcMessage("Parsing Physical Space Basis with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a ReferenceSpaceBasis with the expected dimension."));

    const auto rs = container->get_object<RefSpaceType>(rs_id);
    Assert (rs != nullptr, ExcNullPtr());

    const auto dm_tag = xml_elem->get_single_element("Domain");
    const auto dm_id = dm_tag->get_attribute<Index>("GetFromIgaObjectId");
    AssertThrow (container->is_object<DomainType> (dm_id),
                 ExcMessage("Parsing Physical Space Basis with IgaObjectId=" +
         to_string(object_id) + " the GetFromIgaObjectId does not "
         "correspond to a Domain with the expected dimension."));

    const auto dm = container->get_object<DomainType>(dm_id);
    Assert (dm != nullptr, ExcNullPtr());

    const auto ps = PhysSpaceType::create(rs, dm, transf);

    container->insert_object<PhysSpaceType>(ps, object_id);
}

IGA_NAMESPACE_CLOSE

#endif // XML_IO
