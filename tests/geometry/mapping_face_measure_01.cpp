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
 *  Test for the spherical mapping class.
 *
 *  author: antolin
 *  date: 2014-03-19
 *
 */

#include "../tests.h"

#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/mapping_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/io/writer.h>

template <int dim>
void test_evaluate()
{
    auto grid = CartesianGrid<dim>::create();
    auto map = CylindricalAnnulus::create(grid, 1.0, 2.0, 0.0, 1.0, 0.0, 1.57079632679490);

    QGauss<dim> quad(3);
    auto elem     = map->begin();
    auto elem_end = map->end();
    ValueFlags flag = ValueFlags::face_w_measure|ValueFlags::face_point;
    elem->init_cache(flag, quad);

    std::array<Real, UnitElement<dim>::faces_per_element> face_area ;
    std::fill(face_area.begin(), face_area.end(), 0.0) ;

    for (; elem != elem_end; ++elem)
    {
        if (elem->is_boundary())
        {
            for (Index face_id = 0; face_id < UnitElement<dim>::faces_per_element; ++face_id)
            {
                if (elem->is_boundary(face_id))
                {
                    elem->fill_face_cache(face_id);
                    auto &w_meas = elem->get_face_w_measures(face_id);
                    for (int q = 0; q < w_meas.size(); ++q)
                        face_area[face_id] += w_meas[q] ;
                }
            }
        }
    }

    out << "Dimension " << dim << endl;
    for (Index face_id = 0; face_id < UnitElement<dim>::faces_per_element; ++face_id)
    {
        out << "Area of face " << face_id << " : " << face_area[face_id] << endl;
    }


}

int main()
{
    out.depth_console(10);

//    test_evaluate<2>();
    test_evaluate<3>();

}
