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
 *  Test for the BallFunction class as a mapping
 *
 *  author: pauletti
 *  date: 2014-10-24
 *
 */

#include "../tests.h"

#include <igatools/base/function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_element.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>


template <int dim>
void mapping_values()
{
    using Function = functions::BallFunction<dim>;

    auto flag = NewValueFlags::point | NewValueFlags::value |
                NewValueFlags::gradient |
                NewValueFlags::hessian |
                NewValueFlags::measure|
                NewValueFlags::w_measure;

    auto quad = QUniform<dim>(3);
    auto grid = CartesianGrid<dim>::create();

    auto F = Function::create(grid, IdentityFunction<dim>::create(grid));


    using Mapping   = Mapping<dim, 0>;
    using ElementIt = typename Mapping::ElementIterator;
    Mapping map(F);
    map.reset(flag, quad);

    ElementIt elem(grid, 0);
    ElementIt end(grid, IteratorState::pass_the_end);

    map.template init_cache<dim>(elem);
    for (; elem != end; ++elem)
    {
        map.template fill_cache<dim>(elem, 0);

        out << "Points:" << endl;
        elem->get_points().print_info(out);
        out << endl;
        out << "Values:" << endl;
        elem->template get_values<0, dim>(0).print_info(out);
        out << endl;
        out << "Gradients:" << endl;
        elem->template get_values<1, dim>(0).print_info(out);
        out << endl;
        out << "Hessians:" << endl;
        elem->template get_values<2, dim>(0).print_info(out);
        out << endl;
        out << "Measure:" << endl;
        elem->template get_measures<dim>(0).print_info(out);
        out << endl;
        out << "weight * measure:" << endl;
        elem->template get_w_measures<dim>(0).print_info(out);
        out << endl;
    }
}


int main()
{
    out.depth_console(10);

    mapping_values<1>();
    mapping_values<2>();
    mapping_values<3>();

    return 0;
}
