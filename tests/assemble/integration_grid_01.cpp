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

/*
 *  Test integration grid
 *
 *  author: matrinelli
 *  date: Ma 26, 2015
 *
 */

#include "../tests.h"


#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/geometry/grid_tools.h>


#include <igatools/base/quadrature_lib.h>

using namespace::iga;






template <int dim>
void do_test()
{
    OUTSTART

    using Grid = CartesianGrid<dim>;
    using InterGridMap = grid_tools::InterGridMap<dim>;
    using GridHandler = typename Grid::ElementHandler;


    const auto grid_2 = Grid::create(3);
    const auto grid_3 = Grid::create(4);

    out.begin_item("Grid 2x2:");
    grid_2->print_info(out);
    out.end_item();

    out.begin_item("Grid 3x3:");
    grid_3->print_info(out);
    out.end_item();

    InterGridMap map_elem_grid_integration_to_elem_grid_2;
    InterGridMap map_elem_grid_integration_to_elem_grid_3;
    auto grid_integration = grid_tools::build_cartesian_grid_union(
                                *grid_2,
                                *grid_3,
                                map_elem_grid_integration_to_elem_grid_2,
                                map_elem_grid_integration_to_elem_grid_3);


    out.begin_item("Grid integration:");
    grid_integration->print_info(out);
    out.end_item();


    auto elem_grid_integration = grid_integration->cbegin();
    auto elem_end = grid_integration->cend();

    using MapIdSetId = std::map<int,std::set<int>>;

    MapIdSetId id_elem_2_id_elem_integration ;
    MapIdSetId id_elem_3_id_elem_integration;

    for (; elem_grid_integration != elem_end ; ++elem_grid_integration)
    {
        out.begin_item("Element integration:");
        elem_grid_integration->print_info(out);

        const auto &elem_grid_2 = map_elem_grid_integration_to_elem_grid_2.at(elem_grid_integration);
        out.begin_item("Element grid 2x2:");
        elem_grid_2->print_info(out);
        out.end_item();

        const auto &elem_grid_3 = map_elem_grid_integration_to_elem_grid_3.at(elem_grid_integration);
        out.begin_item("Element grid 3x3:");
        elem_grid_3->print_info(out);
        out.end_item();


        const auto id_elem_integration = elem_grid_integration->get_flat_index();
        const auto id_elem_2 = elem_grid_2->get_flat_index();
        const auto id_elem_3 = elem_grid_3->get_flat_index();

        id_elem_2_id_elem_integration[id_elem_2].insert(id_elem_integration);
        id_elem_3_id_elem_integration[id_elem_3].insert(id_elem_integration);

        out.end_item();
    }


    QGauss<dim> quad_integration(2);
    auto flag = ValueFlags::length;
    auto grid_cache = GridHandler::create(grid_integration);
    grid_cache->reset(flag, quad_integration);


    auto elem_integration = grid_integration->begin();
    grid_cache->init_element_cache(elem_integration);

    out.begin_item("Elements grid 2x2");
    for (const auto &data_elem_2 : id_elem_2_id_elem_integration)
    {
        out << "ID: " << data_elem_2.first <<endl;

        out.begin_item("Elements grid integration");
        for (const auto id_elem_integration : data_elem_2.second)
        {
            out << "elem_integration_id: " << id_elem_integration << endl;
            elem_integration->move_to(id_elem_integration);
            grid_cache->fill_element_cache(elem_integration);

            out << "Vertex 0: " << elem_integration->vertex(0) << endl;

            out << "Lenghts: " << elem_integration->template get_coordinate_lengths<dim>(0) << endl;

        }
        out.end_item();
    }
    out.end_item();


    OUTEND
}








int main()
{
    do_test<2>();

    return 0;
}
