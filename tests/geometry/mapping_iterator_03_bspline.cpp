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
 *  Test for the ig mapping class iterator, geometrical quantities
 *  author: antolin
 *  date: 2014-04-23
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/geometry/mapping_element_accessor.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/io/reader.h>

template <int dim>
void run_test(std::string &file_name)
{
    out << "========== Test (Dimension: " << dim << ") --- begin ========== " << endl;

    // Reading input file.
    auto map = dynamic_pointer_cast<IgMapping<BSplineSpace<dim,dim,1>>>(get_mapping_from_file<dim,0>(file_name));
    map->print_info(out);
    out << endl;

    QTrapez<dim> quad(Real(0.0));
    const auto quad_pts = quad.get_points().get_flat_cartesian_product();
    out << "Quad pts.= " << quad_pts << endl;

    const auto ref_space = map->get_iga_space();
    out << endl;

    //------------------------------------------------------
    out << "Loop using the BSplineElementAccessor" << endl;
    auto sp_elem     = ref_space->begin();
    auto sp_elem_end = ref_space->end();

    sp_elem->init_values(ValueFlags::value, quad);
    for (; sp_elem != sp_elem_end; ++sp_elem)
    {
        sp_elem->fill_values();

        out << "Element id: " << sp_elem->get_flat_index() << endl;

        const auto values = sp_elem->get_basis_values();
        out << "Values = ";
        values.print_info(out);
        out<< endl;
    }
    out << endl;
    //------------------------------------------------------



    //------------------------------------------------------
    out << "Loop using the MappingElementAccessor" << endl;
    auto map_elem     = map->begin();
    auto map_elem_end = map->end();

    map_elem->init_values(ValueFlags::point, quad);
    for (; map_elem != map_elem_end; ++map_elem)
    {
        map_elem->fill_values();
        out << "Element id: " << map_elem->get_flat_index() << endl;

        auto points = map_elem->get_map_values();
        int qp = 0;
        for (auto p : points)
            out << "    Point " << ++qp << ": " << p << endl;
    }
    //------------------------------------------------------
    out << endl;
    out << "========== Test (Dimension: " << dim << ") --- end ========== " << endl;

    out << endl << endl << endl;
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
