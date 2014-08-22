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
 *  Test the IgMapping class on Bspline space
 *  The map is the identity of degree one.
 *  author: pauletti
 *  date: 2013-10-04
 *
 */

#include "../tests.h"

#include <igatools/geometry/ig_mapping.h>
#include <igatools/geometry/mapping_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/io/writer.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>


template <int dim>
void run_test()
{
    typedef BSplineSpace<dim, dim> Space_t;

    const int p = 1;
    auto knots = CartesianGrid<dim>::create(2);
    auto bspline_space = Space_t::create(p, knots);

    vector<Real> control_pts(bspline_space->get_num_basis());
    if (dim == 1)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;
    }
    else if (dim == 2)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
    }
    else if (dim == 3)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

    }


    QGauss<dim> quad(3);

    auto map = IgMapping<Space_t>::create(bspline_space, control_pts);
    ValueFlags flag = ValueFlags::point|ValueFlags::map_gradient;

    auto elem = map->begin();
    elem->init_cache(flag, quad);
    elem->fill_cache();

    auto values = elem->get_map_values();
    auto gradients = elem->get_map_gradients();

    out << "Dim: " << dim << endl;
    out << "Degree: " << p << endl;
    out << "Points: " << endl;
    quad.get_points().get_flat_cartesian_product().print_info(out);
    out << "Values (x1,x2,...):" << endl;
    values.print_info(out);
    out << endl;
    out << "Gradients:" << endl;
    gradients.print_info(out);
    out << endl;


//    //QUniform<dim> quad1(2);
//    QGauss<dim> quad1(2);
//    auto elem1 = map->begin();
//    elem1->init_values(flag, quad1);
//    elem1->fill_values();
//    //elem1->fill_values();
//
////    auto values1 = elem1->get_values();
////    values1.print_info(out);
////    out << endl;
//
//    QUniform<dim> quad2(3);
//    auto elem2 = map->begin();
//    elem2->init_values(flag, quad2);
//    elem2->fill_values();
//    auto values2 = elem2->get_values();
//    values2.print_info(out);
//    out << endl;

    string filename = "bspline_map-" + to_string(dim) + "d";

    Writer<dim> writer(map, 4);
    writer.save(filename);

}

int main()
{
    out.depth_console(10);

    run_test<1>();
    run_test<2>();
    run_test<3>();

}
