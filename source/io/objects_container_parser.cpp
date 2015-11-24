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
    this->parse_grids (xml_elem, container);
    this->parse_spline_spaces (xml_elem, container);
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


IGA_NAMESPACE_CLOSE

#endif // XML_IO
