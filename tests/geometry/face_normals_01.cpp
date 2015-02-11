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
 *  Test of faces normals using a cylindrical annulus mapping.
 *  author: antolin
 *  date: 2014-03-18
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/geometry/push_forward_element_accessor.h>
#include <igatools/io/writer.h>

int main()
{
    out.depth_console(20);


    auto map = CylindricalAnnulus::create(1, 2, 0, 2.0, 0.0, numbers::PI / 4.0);

    const int num_pts = 1 ;
    QGauss<3> quad(num_pts) ;

    auto elem = map->begin();
    ValueFlags flag = ValueFlags::face_normal;
    elem->init_cache(flag, quad);

    for (Index face_id = 0 ; face_id < UnitElement<3>::faces_per_element ; face_id++)
    {
        elem->fill_face_cache(face_id);
        auto normals = elem->get_face_normals(face_id);
        out << "Face: " << face_id << endl;
        out << "  Normal vector: " << endl ;;
        normals.print_info(out);
        out << endl;
    }

    return (0) ;
}

