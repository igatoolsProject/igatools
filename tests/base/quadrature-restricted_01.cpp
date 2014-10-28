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

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/unit_element.h>


template<int dim, int k = dim-1>
void quad_collapse(const int num_pts = 3)
{
    OUTSTART

    QGauss<dim> quad(num_pts);
    out << "Original quadrature:" << endl;
    quad.print_info(out);

    for (auto &id : UnitElement<dim>::template elems_ids<k>())
    {
        auto collapsed_quad = quad.template collapse_to_sub_element<k>(id);
        out << "Quad collapse to subelement: "<< id << endl;
        collapsed_quad.print_info(out);
        out << endl;
    }

    OUTEND
}

int main()
{

    quad_collapse<1>();
    quad_collapse<2>();
    quad_collapse<3>();

    return 0;
}

