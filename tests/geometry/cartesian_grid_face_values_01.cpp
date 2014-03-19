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
 *  Test for the CartesianGrid element iterator
 *  when getting face related values.
 *
 *  author: pauletti
 *  date: Oct 8, 2013
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>

template <int dim>
void run_test1()
{
    out << "========================================================================" << endl;
    out << "Reference faces values <" << dim << ">()" << endl;
    out << endl;

    const int n_knots = 3;
    auto grid = CartesianGrid<dim>::create(n_knots);

    auto elem = grid->begin();
    ValueFlags flag = ValueFlags::face_measure|
                      ValueFlags::face_w_measure|
                      ValueFlags::face_point|
                      ValueFlags::point|
                      ValueFlags::measure|
                      ValueFlags::w_measure;
    elem->init_values(flag, QUniform<dim>(2));
    for (; elem != grid->end(); ++elem)
    {
        out << "Element: "<< elem->get_flat_index() << endl;
        out.push("  ");
        for (int face_id=0; face_id<UnitElement<dim>::faces_per_element; ++face_id)
        {
            elem->fill_face_values(face_id);
            out << "face: " << face_id << endl;
            out.push("  ");
            out << "meas: "<< elem->get_measure(FaceTopology(face_id)) << endl;
            out << "w_meas: "<< endl;
            elem->get_w_measures(FaceTopology(face_id)).print_info(out);
            out << endl;
            out << "points: " << elem->get_points(FaceTopology(face_id)) << endl;
            out.pop();
        }
        out.pop();
    }

    out << "========================================================================" << endl ;
    out << endl ;
}



template <int dim>
void run_test2()
{
    out << "========================================================================" << endl;
    out << "Reference faces values <" << dim << ">()" << endl;
    out << endl;

    TensorSize<dim> n_knots;
    int j=1;
    for (int i = 0 ; i < dim ; ++i)
        n_knots(i) = ++j;

    auto grid = CartesianGrid<dim>::create(n_knots);

    auto elem = grid->begin();
    ValueFlags flag = ValueFlags::face_measure|
                      ValueFlags::face_w_measure|
                      ValueFlags::face_point|
                      ValueFlags::point|
                      ValueFlags::measure|
                      ValueFlags::w_measure;
    elem->init_values(flag, QUniform<dim>(2));
    for (; elem != grid->end(); ++elem)
    {
        out << "Element: "<< elem->get_flat_index() << endl;
        out.push("  ");
        for (int face_id=0; face_id<UnitElement<dim>::faces_per_element; ++face_id)
        {
            elem->fill_face_values(face_id);
            out << "face: " << face_id << endl;
            out.push("  ");
            out << "meas: "<< elem->get_measure(FaceTopology(face_id)) << endl;
            out << "w_meas: "<< endl;
            elem->get_w_measures(FaceTopology(face_id)).print_info(out);
            out << endl;
            out << "points: " << elem->get_points(FaceTopology(face_id)) << endl;
            out.pop();
        }
        out.pop();
    }

    out << "========================================================================" << endl ;
    out << endl ;
}


int main()
{
    out.depth_console(10);

    run_test1<1>();
    run_test1<2>();
    run_test1<3>();

    run_test2<1>();
    run_test2<2>();
    run_test2<3>();

    return  0;
}
