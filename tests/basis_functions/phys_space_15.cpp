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
 *  Test for physical space element values
 *
 *  author: pauletti
 *  date: 2014-11-08
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/identity_function.h>

#include <igatools/basis_functions/space_element_handler.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/geometry/push_forward_element.h>


template<int dim, int codim>
using MapFunc = Function<dim, 0, dim+codim>;


template<int dim, int codim=0>
auto
create_function(shared_ptr<CartesianGrid<dim>> grid)
{

    using Function = functions::LinearFunction<dim, 0, dim+codim>;
    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<dim+codim; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    return Function::create(grid, IdentityFunction<dim>::create(grid), A, b);
}


template <int dim, int order = 0, int range=1, int rank=1, int codim = 0>
void elem_values(const int n_knots = 5, const int deg=1)
{
    OUTSTART
    const int k = dim;
    using RefSpace = BSplineSpace<dim, range, rank>;
    using Space = PhysicalSpace<RefSpace, codim, Transformation::h_grad>;
    using ElementHandler = typename Space::ElementHandler;

    auto grid  = CartesianGrid<dim>::create(n_knots);

    auto ref_space = RefSpace::create(deg, grid);
    auto map_func = create_function(grid);

    auto space = Space::create(ref_space, map_func);

    auto flag = NewValueFlags::none;
    switch (order)
    {
        case 0:
            flag |= NewValueFlags::value;
            break;
        case 1:
            flag |= NewValueFlags::gradient;
            break;
        case 2:
            flag |= NewValueFlags::hessian;
            break;
    }

    auto quad = QGauss<k>(2);

    ElementHandler sp_values(space);
    sp_values.template reset<k> (flag, quad);

    auto elem = space->begin();
    auto end = space->end();
    sp_values.template init_cache<k>(elem);
    for (; elem != end; ++elem)
    {
        sp_values.template fill_cache<k>(elem,0);
        elem->template get_values<0, k>().print_info(out);
    }

    OUTEND
}



int main()
{
    out.depth_console(10);


    elem_values<1>();

    return  0;
}
