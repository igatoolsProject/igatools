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

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
// [quad include]
#include <igatools/base/quadrature_lib.h>
// [quad include]
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

template <int dim>
void loop_on_grid_with_cache()
{
    // [loop as before]
    out << "Traversing the elements of a " << dim << "-dimensional grid." << endl;
    const int n_knots = 3;
    auto grid = CartesianGrid<dim>::create(n_knots);

    using ElementHandler = typename CartesianGrid<dim>::ElementHandler;
    auto elem_handler = ElementHandler::create(grid);

    auto quad = QGauss<dim>(2);
    auto flag = ValueFlags::w_measure;

    elem_handler->template reset<dim>(flag, quad);

    auto elem = grid->begin();
    const auto elem_end = grid->end();
    // [loop as before]
    // [init cache]
    elem_handler->init_element_cache(elem);
    // [init cache]

    for (; elem != elem_end; ++elem)
    {
        // [fill cache]
        elem_handler->fill_element_cache(elem);
        // [fill cache]
        out << "The tensor index of element: " << elem->get_flat_index();
        out << " is: "<< elem->get_tensor_index() << endl;

        // [get meas]
        auto w_meas = elem->template get_w_measures<dim>(0);
        out << "The weighted measure is: ";
        w_meas.print_info(out);
        // [get meas]
        out << endl;
    }
    out << endl;
}


template <int dim>
void loop_on_space_with_cache()
{
    out << "Traversing the elements of a " << dim;
    out << "-dimensional B-spline space." << endl;
    const int n_knots = 3;
    auto grid = CartesianGrid<dim>::create(n_knots);
    const int degree = 2;
    auto space = BSplineSpace<dim>::create(degree, grid);

    using ElementHandler = typename BSplineSpace<dim>::ElementHandler;
    auto elem_handler = ElementHandler::create(space);
    auto quad = QGauss<dim>(1);
    auto flag = ValueFlags::value;

    elem_handler->reset(flag, quad);

    auto elem = space->begin();
    const auto elem_end = space->end();
    elem_handler->init_element_cache(elem);

    for (; elem != elem_end; ++elem)
    {
        elem_handler->fill_element_cache(elem);
        out << "Element: " << elem->get_flat_index();
        out << " has global basis: ";
        elem->get_local_to_global().print_info(out);
        out << endl;
        elem->template get_values<0, dim>(0).print_info(out);
        out<< endl;
    }
    out << endl;
}


int main()
{

    loop_on_grid_with_cache<1>();
    loop_on_grid_with_cache<2>();
    loop_on_grid_with_cache<3>();

    loop_on_space_with_cache<1>();
    loop_on_space_with_cache<2>();
    loop_on_space_with_cache<3>();

    return 0;
}



