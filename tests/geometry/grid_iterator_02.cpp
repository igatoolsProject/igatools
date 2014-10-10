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
 *  using its cache features
 *
 *  author: pauletti
 *  date: Aug 21, 2013
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>


template <int dim>
void run_test()
{
    out << "========================================================================" << endl ;
    out << "Output for function run_test<" << dim << ">()" << endl ;
    out << endl ;


    const int n_knots = 5;
    auto grid = CartesianGrid<dim>::create(n_knots);

    grid->print_info(out);

    out << endl ;


    QGauss<dim> quad(2);
    GridUniformQuadCache<dim> cache(grid, ValueFlags::measure | ValueFlags::w_measure, quad);

    auto elem = grid->begin();
    cache.init_element_cache(elem);
    for (; elem != grid->end(); ++elem)
    {
        out << elem->get_flat_index() << "   ";
        cache.fill_element_cache(elem);
        out << elem->get_measure() << endl;
        elem->get_w_measures().print_info(out);
        out << endl;
    }


    out << "========================================================================" << endl ;
    out << endl ;
}

template <int dim>
void run_test2()
{
    out << "========================================================================" << endl ;
    out << "Output for function run_test2<" << dim << ">()" << endl ;
    out << endl ;


    const int n_knots = 5;
    auto grid = CartesianGrid<dim>::create(n_knots);

    grid->print_info(out);

    out << endl ;

    GridUniformQuadCache<dim> cache(grid, ValueFlags::point, QGauss<dim>(2));

    auto elem = grid->begin();
    cache.init_element_cache(elem);
    for (; elem != grid->end(); ++elem)
    {
        out << elem->get_flat_index() << "   ";
        cache.fill_element_cache(elem);
        elem->get_points().print_info(out);
        out << endl;
    }


    out << "========================================================================" << endl ;
    out << endl ;
}

int main()
{
    out.depth_console(10);
    run_test<1>();
    run_test<2>();
    run_test<3>();

    run_test2<1>();
    run_test2<2>();
    run_test2<3>();

    return  0;
}
