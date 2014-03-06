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

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
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

    auto elem = grid->begin();
    const auto elem_end = grid->end();
    // [loop as before]
    // [init cache]
    QGauss<dim> quad(2);
    ValueFlags fill_flag = ValueFlags::w_measure;
    elem->init_values(fill_flag, quad);
    // [init cache]

    for (; elem != elem_end; ++elem)
    {
        // [fill cache]
        elem->fill_values();
        // [fill cache]
        out << "The center of element: " << elem->get_flat_index();
        out << " is: "<< elem->center() << endl;

        // [get meas]
        auto w_meas = elem->get_w_measures();
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
    auto space = BSplineSpace<dim>::create(grid, degree);

    auto elem = space->begin();
    const auto elem_end = space->end();
    elem->init_values(ValueFlags::value, QGauss<dim>(1));
    for (; elem != elem_end; ++elem)
    {
        elem->fill_values();
        out << "Element: " << elem->get_flat_index();
        out << " has global basis: " << elem->get_local_to_global() << endl;
        elem->get_basis_values().print_info(out);
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



