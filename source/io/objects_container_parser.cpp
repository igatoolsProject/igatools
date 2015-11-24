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

        using ValidGridPtrs = typename ObjectsContainer::ValidGrids;
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
    }
}



template <int dim>
void
ObjectsContainerParser::
parse_grid(const shared_ptr<XMLElement> xml_elem,
           const std::shared_ptr<ObjectsContainer> container) const
{
    Assert (xml_elem->get_name() == "Grid",
            ExcMessage("No Grid xml tag."));

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


IGA_NAMESPACE_CLOSE

#endif // XML_IO
