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
 *  Testing quadrature extension from dim to dim+1
 *  author: pauletti
 *  date: 2014-26 -02
 *
 */

#include "../tests.h"
#include <igatools/base/quadrature.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/unit_element.h>

template <int dim>
void run_test()
{
    const int num_pts = 3;
    out << "======  Extension from "<< dim <<"D to "<<dim+1<<"D  ======" << endl;
    QGauss<dim> quad_surf(num_pts);
    out << "Original quadrature" << endl;
    quad_surf.print_info(out);
    out << endl;
    for (int i = 0; i < UnitElement<dim+1>::faces_per_element; ++i)
    {
        Quadrature<dim+1> quad_extended = quad_surf.get_extension(i);
        quad_surf.print_info(out);
    }
}


int main()
{

    run_test<2>();
    run_test<1>();

    return 0;
}

