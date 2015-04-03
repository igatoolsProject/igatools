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
 *  Test for the push foward using a mapping
 *  with with a linear function
 *  Getting Mapping quantities
 *
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/geometry/push_forward.h>
#include <igatools/geometry/push_forward_element.h>
#include <igatools/geometry/mapping_element.h>

#include <igatools/base/function_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_element.h>

template<int dim, int codim>
void test()
{
    const int space_dim = dim + codim;
    using Function = functions::LinearFunction<dim, codim, space_dim>;
    using PForward  = PushForward<Transformation::h_grad, dim, codim>;

    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<space_dim; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    auto quad = QGauss<dim>(2);
    auto grid = CartesianGrid<dim>::create(3);

    auto flag = ValueFlags::point | ValueFlags::value
                | ValueFlags::gradient | ValueFlags::hessian|
                ValueFlags::measure |
                ValueFlags::w_measure;
    auto F = Function::create(grid, IdentityFunction<dim>::create(grid), A, b);
    PForward pf(F);

    pf.template reset<dim>(flag, quad);

    auto elem = pf.begin();
    auto end  = pf.end();

    pf.template init_cache<dim>(*elem);

    for (; elem != end; ++elem)
    {
        pf.template fill_cache<dim>(*elem, 0);
        elem->get_points().print_info(out);
        out << endl;
        elem->template get_values<_Value,dim>(0).print_info(out);
        out << endl;
        elem->template get_values<_Gradient,dim>(0).print_info(out);
        out << endl;
        elem->template get_values<_Hessian,dim>(0).print_info(out);
        out << endl;
        elem->template get_measures<dim>(0).print_info(out);
        out << endl;
        elem->template get_w_measures<dim>(0).print_info(out);
        out << endl;
    }


}


int main()
{
    test<2,0>();
//    test<3,3>();

    return 0;
}

