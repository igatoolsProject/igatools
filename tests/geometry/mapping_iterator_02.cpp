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
 *  Test for the linear mapping class iterator, geometrical quantities
 *  author: pauletti
 *  date: 2013-09-26
 *
 */

#include "../tests.h"

#include <igatools/geometry/mapping_lib.h>
#include <igatools/base/quadrature_lib.h>

#include <igatools/geometry/mapping_element_accessor.h>

template <int dim, int codim>
void test_iterator()
{
    const int space_dim = dim + codim;
    auto grid = CartesianGrid<dim>::create(3);
    Derivatives<dim,space_dim,1,1> A;
    Point<space_dim> b;
    //Dilation
    for (int i=0; i<dim; ++i)
        A[i][i] = 3.;
    //Traslation
    for (int i=0; i<dim; ++i)
        b[i] = 1;

    auto map = LinearMapping<dim,codim>::create(grid,A,b);

    out << "Linear mapping" << "<" << dim << "," << space_dim << ">" << endl;
    out << "A =" << endl << A << endl;
    out << "b =" << b << endl << endl;

    QTrapez<dim> quad;

    auto elem = map->begin();

    ValueFlags flag = ValueFlags::w_measure;
    flag |= ValueFlags::measure| ValueFlags::map_gradient;
    elem->init_values(flag, quad);

    elem->fill_values();

//   auto values = elem->get_normals();
//    auto dets = elem->get_dets_map();
    auto wdets = elem->get_w_measures();
//    auto gradients = elem->get_gradients_map();
//    auto hessians = elem->get_hessians_map();
//
//    out << "x = " << endl << quad.get_points().get_flat_cartesian_product() << endl;
//    out << "F(x)     = " << endl;
//    values.print_info(out);
//    out << endl;
//    out << "dets    = " << endl;
//    dets.print_info(out);
//    out << endl;
    out << "wdets    = " << endl;
    wdets.print_info(out);
    out << endl;
//    out << "Dx    = " << endl;
//    gradients.print_info(out);
//    out << endl;
//    out << "D2x   = " << endl;
//    hessians.print_info(out);
    out << endl;
}

int main()
{
    out.depth_console(10);
    test_iterator<1,0>();
//   test_iterator<2,2>();
//   test_iterator<3,3>();
//    test_evaluate<2,3>();
//    test_evaluate<1,2>();
//    test_evaluate<1,3>();

}
