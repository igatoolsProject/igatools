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
 *  Test for LinearFunction
 *  author: pauletti
 *  date: Oct 12, 2014
 */

#include "../tests.h"

#include <igatools/base/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/function_element.h>



template<int dim, int codim, int range>
void test(shared_ptr<Function<dim, codim, range>> F)
{
    using ElementIterator = typename  Function<dim, codim, range>::ElementIterator;
    ElementIterator elem = F->begin();
    ElementIterator end = F->end();

    F->init_cache(elem, Int<dim>());
    for (; elem != end; ++elem)
    {
        F->fill_cache(elem, Int<dim>(),0);
        elem->get_points().print_info(out);
        out << endl;
        elem->template get_values<0, dim>(0).print_info(out);
        out << endl;
        elem->template get_values<1, dim>(0).print_info(out);
        out << endl;
        elem->template get_values<2, dim>(0).print_info(out);
        out << endl;
    }
}


template<int dim, int codim, int range>
void create_fun()
{
    using Function = functions::LinearFunction<dim, codim, range>;

    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<range; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    auto flag = ValueFlags::point | ValueFlags::value | ValueFlags::gradient |
                ValueFlags::hessian;
    auto quad = QGauss<dim>(2);
    auto grid = CartesianGrid<dim>::create(3);
    auto F = Function::create(grid, IdentityFunction<dim>::create(grid), A, b);
    F->reset(flag, quad);
    test<dim, codim, range>(F);
}




int main()
{
    create_fun<2, 0, 2>();

    return 0;
}

