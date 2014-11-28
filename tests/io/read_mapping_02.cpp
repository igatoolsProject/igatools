//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
 *  Test for the ig mapping class iterator, geometrical quantities using a NURBS mapping
 *  author: antolin
 *  date: 2014-04-23
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/ig_function.h>
#include <igatools/basis_functions/nurbs_element.h>
#include <igatools/io/reader.h>

template <int dim>
void run_test(std::string &file_name)
{
    OUTSTART

    // Reading input file.
    using RefSpace = NURBSSpace<dim,dim,1>;
    auto map = dynamic_pointer_cast<IgFunction<RefSpace> >(get_mapping_from_file<dim,0>(file_name));
    out.begin_item("IgFunction infos:");
    map->print_info(out);
    out << endl;

    QTrapez<dim> quad(Real(0.0));
    const auto quad_pts = quad.get_points().get_flat_cartesian_product();
    out.begin_item("Quad pts.:");
    quad_pts.print_info(out);
    out.end_item();
    out << endl;


    const auto ref_space = map->get_iga_space();

    //------------------------------------------------------
    out.begin_item("Loop using the NURBSElement");
    using ElemHandler = typename RefSpace::ElementHandler;
    ElemHandler sp_elem_handler(ref_space);
    sp_elem_handler.reset(NewValueFlags::value,quad);

    auto sp_elem     = ref_space->begin();
    auto sp_elem_end = ref_space->end();

    sp_elem_handler.template init_cache<dim>(sp_elem);
    for (; sp_elem != sp_elem_end; ++sp_elem)
    {
        sp_elem_handler.template fill_cache<dim>(sp_elem,0);

        out << "Element id: " << sp_elem->get_flat_index() << endl;

        const auto &values = sp_elem->template get_values<0,dim>(0);
        out << "Values = ";
        values.print_info(out);
        out<< endl;
    }
    out.end_item();
    out << endl;
    //------------------------------------------------------


    //------------------------------------------------------
    out.begin_item("Loop using the FunctionElement");
    map->reset(NewValueFlags::point, quad);

    auto map_elem     = map->begin();
    auto map_elem_end = map->end();

    const auto topology = Int<dim>();

    map->init_cache(*map_elem,topology);

    for (; map_elem != map_elem_end; ++map_elem)
    {
        map->fill_cache(*map_elem,0,topology);
        out << "Element id: " << map_elem->get_flat_index() << endl;

        const auto &points = map_elem->get_points();
        int qp = 0;
        for (auto p : points)
            out << "    Point " << ++qp << ": " << p << endl;
    }
    out.end_item();
    out << endl;
    //------------------------------------------------------
    OUTEND
}


int main()
{

    std::string file_name = "cube_1D.xml";
    run_test<1>(file_name);

    file_name = "cube_2D.xml";
    run_test<2>(file_name);
    file_name = "cube_3D.xml";
    run_test<3>(file_name);
//*/
    return (0) ;
}
