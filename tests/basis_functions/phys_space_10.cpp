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
 *  Test for "New" physical space
 *
 *  author: pauletti
 *  date: Oct 08, 2014
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


template <int dim, int range=1, int rank=1, int codim = 0>
void cache_init(const NewValueFlags flag,
                const int n_knots = 5, const int deg=1)
{
    OUTSTART

    using RefSpace = NewBSplineSpace<dim, range, rank>;
    using Space    = PhysicalSpace<RefSpace, codim, Transformation::h_grad>;
    using ElementHandler = typename Space::ElementHandler;
    auto grid      = CartesianGrid<dim>::create(n_knots);
    auto ref_space = RefSpace::create(deg, grid);

    using Function = functions::LinearFunction<dim, 0, dim+codim>;
    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<Space::space_dim; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    auto quad = QGauss<dim>(2);
    auto map_func = Function::create(grid, IdentityFunction<dim>::create(grid), A, b);
    auto space = Space::create(ref_space, map_func);


    ElementHandler sp_values(space);
    sp_values.template reset<dim> (flag, quad);
    sp_values.print_info(out);

    OUTEND
}



template <int dim, int range=1, int rank=1, int codim = 0>
void cache_init_elem(const NewValueFlags flag,
                     const int n_knots = 5, const int deg=1)
{
    const int k = dim;
    OUTSTART

    using RefSpace = NewBSplineSpace<dim, range, rank>;
    using Space = PhysicalSpace<RefSpace, codim, Transformation::h_grad>;
    using ElementHandler = typename Space::ElementHandler;

    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto ref_space = RefSpace::create(deg, grid);

    using Function = functions::LinearFunction<dim, 0, dim+codim>;
    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<Space::space_dim; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    auto quad = QGauss<dim>(2);
    auto map_func = Function::create(grid, IdentityFunction<dim>::create(grid), A, b);
    auto space = Space::create(ref_space, map_func);

    ElementHandler sp_values(space);
    sp_values.template reset<dim> (flag, quad);

    auto elem = space->begin();
    sp_values.template init_cache<k>(elem);
    elem->print_cache_info(out);

    OUTEND
}


template <int dim, int range=1, int rank=1, int codim = 0>
void cache_fill_elem(const NewValueFlags flag,
                     const int n_knots = 5, const int deg=1)
{
    OUTSTART

    const int k = dim;
    using RefSpace = NewBSplineSpace<dim, range, rank>;
    using Space = PhysicalSpace<RefSpace, codim, Transformation::h_grad>;
    using ElementHandler = typename Space::ElementHandler;

    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto ref_space = RefSpace::create(deg, grid);

    using Function = functions::LinearFunction<dim, 0, dim+codim>;
    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<Space::space_dim; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    auto quad = QGauss<dim>(2);
    auto map_func = Function::create(grid,IdentityFunction<dim>::create(grid), A, b);
    auto space = Space::create(ref_space, map_func);

    ElementHandler sp_values(space);
    sp_values.template reset<dim> (flag, quad);

    auto elem = space->begin();
    auto end = space->end();
    sp_values.template init_cache<k>(elem);
    for (; elem != end; ++elem)
    {
        sp_values.template fill_cache<k>(elem, 0);
        elem->print_cache_info(out);
    }

    OUTEND
}



template <int dim, int range=1, int rank=1, int codim = 0>
void cache_get_elem_values(const NewValueFlags flag,
                           const int n_knots = 5, const int deg=1)
{
    OUTSTART
    const int k = dim;
    using RefSpace = NewBSplineSpace<dim, range, rank>;
    using Space = PhysicalSpace<RefSpace, codim, Transformation::h_grad>;
    using ElementHandler = typename Space::ElementHandler;

    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto ref_space = RefSpace::create(deg, grid);

    using Function = functions::LinearFunction<dim, 0, dim+codim>;
    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<Space::space_dim; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    auto quad = QGauss<dim>(2);
    auto map_func = Function::create(grid, IdentityFunction<dim>::create(grid), A, b);
    auto space = Space::create(ref_space, map_func);

    ElementHandler sp_values(space);
    sp_values.template reset<dim> (flag, quad);

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

    cache_init<1>(NewValueFlags::value);
    cache_init_elem<1>(NewValueFlags::value);
    cache_fill_elem<1>(NewValueFlags::value);
    cache_get_elem_values<1>(NewValueFlags::value);

    return  0;
}
