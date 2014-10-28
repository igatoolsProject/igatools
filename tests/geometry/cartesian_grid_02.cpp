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
 *  Test for CartesianGrid face boundary ids (get and set)
 *
 *  author: pauletti
 *  date: 2012-11-27
 *  QA: v0.2 (pauletti, 2013-10-25)
 */

#include "../tests.h"
#include <igatools/geometry/cartesian_grid.h>

template <int dim>
void test_face_id()
{
    auto cartesian_grid = CartesianGrid<dim>::create();
    out << endl << "CartesianGrid Dimension: " << dim << endl;

    for (int j = 0; j < UnitElement<dim>::n_faces; ++j)
    {
        out << "Face number: " << j << endl;
        out << "Face boundary id: " << cartesian_grid->get_boundary_id(j) << endl;
    }

    for (int j = 0; j < UnitElement<dim>::n_faces; ++j)
        cartesian_grid->set_boundary_id(j,j);

    for (int j = 0; j < UnitElement<dim>::n_faces; ++j)
    {
        out << "Face number: " << j << endl;
        out << "Face boundary id: " << cartesian_grid->get_boundary_id(j) << endl;
    }
}



int main()
{
//  test_face_id<0>();
    test_face_id<1>();
    test_face_id<2>();
    test_face_id<3>();

    return 0;
}
