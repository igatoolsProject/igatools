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
 *  Testing quadrature restriction
 *  author: pauletti
 *  date: 2013-04-02
 *
 */

#include "../tests.h"
#include <igatools/base/quadrature.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/unit_element.h>


template<int dim>
void do_test()
{
    out<< "Dimension: " << dim << endl;
    const int num_pts = 3;

    QGauss<dim> quad(num_pts);
    out << "Original quadrature:" << endl;
    quad.print_info(out);

    for (int face_id = 0; face_id<UnitElement<dim>::faces_per_element; ++face_id)
    {
        auto new_quad = restricted_quad(quad, face_id);
        out << "Restricted quadrature to face: "<< face_id << endl;
        new_quad.print_info(out);
    }

}

int main()
{

    do_test<1>();
    do_test<2>();
    do_test<3>();

    return 0;
}

