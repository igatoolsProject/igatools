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
 *  Test for the spherical mapping class.
 *
 *  author: pauletti
 *  date: 2013-02-20
 *
 */

#include "../tests.h"

#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/mapping_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/io/writer.h>

template <int dim>
void test_evaluate()
{
    auto grid = CartesianGrid<dim>::create();
    auto map = BallMapping<dim>::create(grid);

    QGauss<dim> quad(1);
    ValueFlags flag = ValueFlags::point|ValueFlags::map_gradient
                      |ValueFlags::map_hessian;

    auto elem = map->begin();
    elem->init_cache(flag, quad);
    elem->fill_cache();

    auto values    = elem->get_map_values();
    auto gradients = elem->get_map_gradients();
    auto hessians  = elem->get_map_hessians();

    out << "Points: (r,phi1,...,phi_n) :" << endl;
    out << quad.get_points().get_flat_cartesian_product() << endl;

    out << "Values (x1,x2,...):" << endl;
    values.print_info(out);

    out << "Gradients (x1,x2,...):" << endl;
    gradients.print_info(out);

    out << "Hessians:" << endl;
    hessians.print_info(out);

    string filename = "spherical_map-" + to_string(dim) + "d";
    Writer<dim> writer(map, 10);
    writer.save(filename);
}


int main()
{
    out.depth_console(10);

    test_evaluate<1>();
    test_evaluate<2>();
    test_evaluate<3>();

    return 0;
}
