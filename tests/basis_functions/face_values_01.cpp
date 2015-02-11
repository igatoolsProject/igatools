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
 *  Test for the PhysicalSpace element iterator to get the face values
 *  author: antolin, pauletti
 *  date: 2014-02-10
 */


#include "../tests.h"

#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/unit_element.h>


template<int dim>
void run_test(const int num_knots, const int p)
{
    using ref_space_t = BSplineSpace<dim>;
    using pf_t = PushForward<Transformation::h_grad,dim>;
    using space_t = PhysicalSpace<ref_space_t, pf_t>;

    auto grid = CartesianGrid<dim>::create(num_knots);
    auto ref_space = ref_space_t::create(p, grid);
    auto map = IdentityMapping<dim>::create(grid);
    auto space = space_t::create(ref_space, pf_t::create(map));

    const int num_faces = UnitElement<dim>::faces_per_element;
    QGauss<dim> quad(2);

    ValueFlags flag = ValueFlags::face_value;

    auto elem = space->begin();
    auto end =  space->end();

    elem->init_cache(flag, quad);

    for (; elem != end; ++elem)
    {
        if (elem->is_boundary())
        {
            out << "Element" << elem->get_flat_index() << endl;

            for (Index face = 0; face < num_faces; ++face)
            {
                if (elem->is_boundary(face))
                {
                    out << "Face " << face << endl;
                    elem->fill_face_cache(face);
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
    out.depth_console(20);

    const int p = 1;
//    const int num_knots = 2;
//    run_test<2>(num_knots, p);

    for (int num_knots = 2; num_knots<4; ++num_knots)
    {
        run_test<2>(num_knots, p);
        run_test<3>(num_knots, p);
    }

    return 0;
}

