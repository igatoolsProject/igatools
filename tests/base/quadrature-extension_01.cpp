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
 *  Testing quadrature extension from dim to dim+1
 *  author: pauletti
 *  date: 2014-26 -02
 *
 */

#include "../tests.h"
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/unit_element.h>

template <int k, int dim = k+1>
void quad_extension(const int num_pts = 2)
{
    OUTSTART

    QGauss<k> sub_quad(num_pts);
    out << "Original quadrature" << endl;
    sub_quad.print_info(out);
    out << endl;

    for (auto &i : UnitElement<dim>::template elems_ids<k>())
    {
        out << "Extended to subelement: " << i << endl;
        Quadrature<dim> vol_quad = extend_sub_elem_quad<k, dim>(sub_quad, i);
        vol_quad.print_info(out);
    }
    OUTEND
}


int main()
{

    quad_extension<0>();
    quad_extension<1>();
    quad_extension<2>();

    return 0;
}

