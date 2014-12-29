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
 *  Test for IgMapping class without the use of the cache.
 *  The output is the same of the test geometry/bspline_mapping_04
 *  author: martinelli
 *  date: 2014-05-21
 *
 */

#include "../tests.h"

#include <igatools/geometry/ig_mapping.h>
#include <igatools/geometry/mapping_element_accessor.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>

#include <igatools/base/quadrature_lib.h>
#include <igatools/io/writer.h>
#include <igatools/linear_algebra/dof_tools.h>

template <int dim>
void test_evaluate()
{

    const int rank = (dim==1) ? 0 : 1 ;

    //----------------------------------------------------------------------------------------------
    out << "Dim: " << dim << endl ;
    int n_knots = 2;
    CartesianProductArray<Real , dim> coord ;
    for (int i = 0; i < dim; ++i)
    {
        vector<Real> tmp_coord;
        for (int j = 0; j < n_knots; ++j)
            tmp_coord.push_back(j);
        coord.copy_data_direction(i,tmp_coord);
    }


    int p = 2 ;


    auto knots = CartesianGrid<dim>::create(coord);

    typedef BSplineSpace<dim,dim,rank> Space_t ;

    shared_ptr< Space_t > bspline_space = Space_t::create(p, knots)  ;
    //----------------------------------------------------------------------------------------------


    //----------------------------------------------------------------------------------------------
    vector<Real> control_pts(bspline_space->get_num_basis()) ;

    if (dim == 2)
    {
        int id = 0 ;

        // x coords
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.5 ;
        control_pts[id++] = 1.5 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 2.0 ;
        control_pts[id++] = 2.0 ;


        // y coords
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.5 ;
        control_pts[id++] = 1.5 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 2.0 ;
        control_pts[id++] = 2.0 ;
        control_pts[id++] = 0.0 ;
    }
    else if (dim == 3)
    {
        int id = 0 ;

        // x coords
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.5 ;
        control_pts[id++] = 1.5 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 2.0 ;
        control_pts[id++] = 2.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.5 ;
        control_pts[id++] = 1.5 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 2.0 ;
        control_pts[id++] = 2.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.5 ;
        control_pts[id++] = 1.5 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 2.0 ;
        control_pts[id++] = 2.0 ;


        // y coords
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.5 ;
        control_pts[id++] = 1.5 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 2.0 ;
        control_pts[id++] = 2.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.5 ;
        control_pts[id++] = 1.5 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 2.0 ;
        control_pts[id++] = 2.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.5 ;
        control_pts[id++] = 1.5 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 2.0 ;
        control_pts[id++] = 2.0 ;
        control_pts[id++] = 0.0 ;

        // z coords
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 0.5 ;
        control_pts[id++] = 0.5 ;
        control_pts[id++] = 0.5 ;
        control_pts[id++] = 0.5 ;
        control_pts[id++] = 0.5 ;
        control_pts[id++] = 0.5 ;
        control_pts[id++] = 0.5 ;
        control_pts[id++] = 0.5 ;
        control_pts[id++] = 0.5 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

    }



    auto map = IgMapping<Space_t>::create(bspline_space, control_pts) ;


    QGauss<dim> quad(3);
    const auto points = quad.get_points().get_flat_cartesian_product();

    auto elem = map->begin();

    auto values = elem->evaluate_values_at_points(points);
    auto gradients = elem->evaluate_gradients_at_points(points);
    auto hessians = elem->evaluate_hessians_at_points(points);

    out << "Values : ";
    values.print_info(out);
    out << endl;

    out << "Gradients : ";
    gradients.print_info(out);
    out << endl;

    out << "Hessians : ";
    hessians.print_info(out);
    out << endl;

    string filename = "bspline_map-" + to_string(dim) + "d";
    Writer<dim> writer(map, 4);
    writer.save(filename) ;

}

int main()
{
    out.depth_console(10);

    test_evaluate<2>();
    test_evaluate<3>();

}
