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
 *  Test for the push forward using a mapping
 *  with with a linear function
 *
 *  Performing transforms
 *
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/geometry/new_push_forward.h>
#include <igatools/geometry/push_forward_element.h>
#include <igatools/geometry/mapping_element.h>

#include <igatools/base/function_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_element.h>

#include <igatools/basis_functions/new_bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>

template<int dim, int codim>
void test()
{
    const int space_dim = dim + codim;
    using Function = functions::LinearFunction<dim, codim, space_dim>;
    using Space = BSplineSpace<dim>;
    using PForward  = NewPushForward<Transformation::h_grad, dim, codim>;

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

    auto space  = Space::create(1, grid);
    auto sp_flag = NewValueFlags::value | NewValueFlags::gradient | NewValueFlags::hessian;
    typename Space::ElementHandler sp_values(space);
    sp_values.template reset<dim>(sp_flag, quad);
    auto sp_elem = space->begin();


    auto F = Function::create(grid, IdentityFunction<dim>::create(grid), A, b);
    PForward pf(F);

    auto flag = NewValueFlags::point |
                NewValueFlags::tran_value |
                NewValueFlags::tran_gradient|
                NewValueFlags::tran_hessian;

    pf.template reset<dim>(flag, quad);

    auto elem = pf.begin();
    auto end  = pf.end();


    pf.template init_cache<dim>(elem.get_accessor());
    sp_values.template init_cache<dim>(sp_elem);

    for (; elem != end; ++elem, ++sp_elem)
    {
        sp_values.template fill_cache<dim>(sp_elem, 0);
        pf.template fill_cache<dim>(elem.get_accessor(), 0);

        const auto &ref_values = sp_elem->template get_values<0,dim>(0);
        ValueTable<typename PForward::template PhysValue<Space::range, Space::rank>>
        values(ref_values.get_num_functions(), ref_values.get_num_points());
        elem->template transform_0<Space::range, Space::rank>
        (ref_values, values);

        ref_values.print_info(out);
        out << endl;
        values.print_info(out);
        out << endl;

        const auto &ref_der_1 = sp_elem->template get_values<1,dim>(0);
        ValueTable<typename PForward::template PhysDerivative<Space::range, Space::rank, 1>>
        gradients(ref_values.get_num_functions(), ref_values.get_num_points());
        elem->template transform_1<Space::range, Space::rank, dim>
        (std::make_tuple(ref_values, ref_der_1), values, gradients, 0);

        ref_der_1.print_info(out);
        out << endl;
        gradients.print_info(out);
        out << endl;


        const auto &ref_der_2 = sp_elem->template get_values<2,dim>(0);
        ValueTable<typename PForward::template PhysDerivative<Space::range, Space::rank, 2>>
        hessians(ref_values.get_num_functions(), ref_values.get_num_points());
        elem->template transform_2<Space::range, Space::rank, dim>
        (std::make_tuple(ref_values, ref_der_1, ref_der_2),
         std::make_tuple(values, gradients),
         hessians, 0);

        ref_der_2.print_info(out);
        out << endl;
        hessians.print_info(out);
        out << endl;

    }


}


int main()
{



    test<2,0>();
//    test<3,3>();

    return 0;
}

