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
 *  Test for the BSplineSpace element iterator to get the face values
 *  author: pauletti
 *  date: 2013-10-23
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/unit_element.h>

template<int dim>
void run_test(const int num_knots, const int p)
{
    using space_t = BSplineSpace<dim>;

    auto grid = CartesianGrid<dim>::create(num_knots);
    auto space = space_t::create(p, grid);

    const int num_faces = UnitElement<dim>::n_faces;
    QGauss<dim> quad(2);

    ValueFlags flag = ValueFlags::face_value|ValueFlags::value;

    auto elem = space->begin();
    auto end =  space->end();
    elem->init_cache(flag, quad);

    for (; elem != end; ++elem)
    {
        if (elem->is_boundary())
        {
            elem->fill_cache();
            out << "Element" << elem->get_flat_index() << endl;
            elem->get_basis_values().print_info(out);

            for (int face = 0; face < num_faces; ++face)
            {
                if (elem->is_boundary(face))
                {
                    elem->fill_face_cache(face);
                    out << "Face " << face << endl;
                    out << "values: " << endl;
                    elem->get_face_basis_values(face).print_info(out);
                    out << endl;
                }
            }
        }
    }
    out << endl;
}




int main()
{
    out.depth_console(10);
    run_test<2>(2,1);
    //do_test<2>(3,1);

    return 0;
}

