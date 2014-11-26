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
 *  Test for the SphericalFunction class as a mapping
 *
 *  author: pauletti
 *  date: 2014-10-24
 *
 */

#include "../tests.h"

#include <igatools/base/function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_element.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/mapping_element.h>


template <int dim>
void principal_curvatures()
{
    OUTSTART

    using Function = functions::SphereFunction<dim>;

    auto flag = NewValueFlags::point |  NewValueFlags::value |
                NewValueFlags::curvature;

    auto quad = QUniform<dim>(3);

    BBox<dim> box;
    for (int i=0; i<dim-1; ++i)
        box[i] = {0.+M_PI/8, M_PI-M_PI/8};
    if (dim>=1)
        box[dim-1] = {0., M_PI};
    auto grid = CartesianGrid<dim>::create(box, 2);

    auto F = Function::create(grid, IdentityFunction<dim>::create(grid));


    using Mapping   = NewMapping<dim, 1>;
    Mapping map(F);
    map.reset(flag, quad);

    auto elem = map.begin();
    auto end = map.end();

    map.template init_cache<dim>(elem);
    for (; elem != end; ++elem)
    {
        map.template fill_cache<dim>(elem, 0);
        out << "Normals:" << endl;
        elem->get_external_normals().print_info(out);
        out << endl;

        out << "Curvature:" << endl;
        elem->get_principal_curvatures().print_info(out);
        out << endl;
    }

    OUTEND
}


int main()
{
    out.depth_console(10);

    principal_curvatures<1>();
    principal_curvatures<2>();

    return 0;
}
