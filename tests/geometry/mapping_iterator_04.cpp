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
 *  Test for the ig mapping class iterator, geometrical quantities
 *  author: antolin
 *  date: 2014-04-26
 *
 */


#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/geometry/mapping_element_accessor.h>
#include <igatools/geometry/ig_mapping.h>

template<int dim>
void run_test()
{
    out << "========== Test (Dimension: " << dim << ") --- begin ========== " << endl;

    const int range = dim;
    const int rank =  1;
    const int degree = 2;
    const int face_id = 1;
    auto knots = CartesianGrid<dim>::create();
    typedef BSplineSpace< dim, range, rank > RefSpace_t;
    auto ref_space = RefSpace_t::create(knots, degree);

    vector<Real> control_pts(ref_space->get_num_basis());
    if (dim == 2)
    {
        control_pts[0] = 0.0;
        control_pts[1] = 0.5;
        control_pts[2] = 1.0;
        control_pts[3] = 0.0;
        control_pts[4] = 0.5;
        control_pts[5] = 1.0;
        control_pts[6] = 0.0;
        control_pts[7] = 0.5;
        control_pts[8] = 1.0;

        control_pts[9] = 0.0;
        control_pts[10] = 0.0;
        control_pts[11] = 0.0;
        control_pts[12] = 0.5;
        control_pts[13] = 0.75;
        control_pts[14] = 1.0;
        control_pts[15] = 1.0;
        control_pts[16] = 1.5;
        control_pts[17] = 2.0;
    }
    else if (dim == 3)
    {
        AssertThrow(false,ExcNotImplemented());
    }
    else if (dim == 1)
    {
        AssertThrow(false,ExcNotImplemented());
    }

    auto map = IgMapping<RefSpace_t>::create(ref_space, control_pts);
    map->refine_h_direction(0, 2);
    map->print_info(out);
    out << endl;

    shared_ptr<Quadrature<dim>> quad = std::make_shared<Quadrature<dim>>(QTrapez<dim>());

    auto elem     = map->begin();
    auto elem_end = map->end();
    ValueFlags fill_flags = ValueFlags::map_value | ValueFlags::face_point;
    elem->init_values(fill_flags, *quad);


    vector<Point<dim>> unit_points_face(2);
    unit_points_face[0][0] = 1.0;
    unit_points_face[0][1] = 0.0;
    unit_points_face[1][0] = 1.0;
    unit_points_face[1][1] = 1.0;
//*/
    const auto unit_points = quad->get_points().get_flat_cartesian_product();
    int point_id = 0;
    for (const auto &pt : unit_points)
        out << "Point[" << point_id++ << "]= " << pt <<endl;

    out << "Loop using the MappingElementAccessor" << endl;

    for (; elem != elem_end; ++elem)
    {
        out << "----" << endl;
        out << "Element id " << elem->get_flat_index() << endl;
        elem->print_info(out);

        out << "Points evaluated at the element using the cache" << endl;
        elem->fill_values();
        auto points = elem->get_map_values();
        points.print_info(out);
        out << endl;

        out << "Points evaluated at the element using evaluate_values_at_points()" << endl;
        auto points_no_cache = elem->evaluate_values_at_points(unit_points);
        points_no_cache.print_info(out);
        out << endl;

        out << "Points evaluated at the face id " << face_id << " using the cache" << endl;
        elem->fill_face_values(face_id);
        auto points_face = elem->get_face_values(face_id);
        points_face.print_info(out);
        out << endl;

        out << "Points evaluate at the face id " << face_id << " using evaluate_values_at_points()" << endl;
        auto points_face_no_cache = elem->evaluate_values_at_points(unit_points_face);
        points_face_no_cache.print_info(out);
        out << "----" << endl;
        out << endl;
        //*/
    }
    out << "========== Test (Dimension: " << dim << ") --- end ========== " << endl;
}


int main()
{
    run_test<2>();
//     run_test<3>();
    return (0);
}
