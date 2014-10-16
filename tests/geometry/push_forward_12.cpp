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
 *  Test for linear mapping class
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/geometry/new_push_forward.h>
#include <igatools/geometry/push_forward_element.h>
#include <igatools/geometry/new_mapping_element_accessor.h>
#include <igatools/../../source/geometry/grid_forward_iterator.cpp>
#include <igatools/base/function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_element.h>

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/basis_functions/bspline_uniform_quad_cache.h>

template<int dim, int codim>
void test()
{
    const int space_dim = dim + codim;
    using Function = functions::LinearFunction<dim, codim, space_dim>;
    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<space_dim; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    auto flag = ValueFlags::point | ValueFlags::value | ValueFlags::gradient |
                ValueFlags::hessian |ValueFlags::measure|ValueFlags::w_measure;
    auto quad = QGauss<dim>(2);
    auto grid = CartesianGrid<dim>::create(3);


    using Space = BSplineSpace<dim>;
    auto space  = Space::create(1, grid);
    typename Space::UniformQuadCache cache(space, ValueFlags::value, quad);
    auto sp_elem = space->begin();
    cache.init_element_cache(sp_elem);

    auto F = Function::create(grid, flag, quad, A, b);
    using PForward  = NewPushForward<Transformation::h_grad, dim, codim>;
    using ElementIt = typename PForward::ElementIterator;
    PForward pf(F, flag, quad);


    ElementIt elem(grid, 0);
    ElementIt end(grid, IteratorState::pass_the_end);


    pf.init_element(elem);
    for (; elem != end; ++elem)
    {
        cache.fill_element_cache(sp_elem);
        pf.fill_element(elem);
        elem->get_points().print_info(out);
        out << endl;
        elem->get_values().print_info(out);
        out << endl;
        elem->get_gradients().print_info(out);
        out << endl;
        elem->get_hessians().print_info(out);
        out << endl;
        elem->get_measures().print_info(out);
        out << endl;
        elem->get_w_measures().print_info(out);
        out << endl;

        const auto &ref_values = sp_elem->get_basis_values();
        ValueTable<typename PForward::template PhysValue<Space::range, Space::rank>>
        values(ref_values.get_num_functions(), ref_values.get_num_points());
        elem->template transform_0<Space::range, Space::rank>
        (ref_values, values);

        ref_values.print_info(out);
        out << endl;
        values.print_info(out);
        out << endl;
    }


}


int main()
{



    test<2,0>();
//    test<3,3>();

    return 0;
}

