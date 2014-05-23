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
#include <igatools/geometry/mapping_element_accessor.h>
#include <igatools/io/reader.h>

template <int dim>
void run_test(std::string& file_name)
{
    out << "Dimension: " << dim << endl;

    // Reading input file.
    auto map = get_mapping_from_file<dim,0>(file_name);

    QTrapez<dim> quad(Real(0.0));

    auto elem     = map->begin();
    auto elem_end = map->end();

    ValueFlags flag = ValueFlags::point;
    elem->init_values(flag, quad);
    for (; elem != elem_end; ++elem)
    {
        out << "Element: " << elem->get_flat_index() << endl;
        elem->fill_values();
        auto points = elem->get_values();
        int qp = 0;
        for(auto p : points)
            out << "    Point " << ++qp << ": " << p << endl;
    }

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

    return (0) ;
}
