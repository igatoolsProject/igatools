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
 *  Test for BSplineSpace element iterator using
 *  the uniform quad global cache
 *  Evaluates values, gradients and derivatives
 *
 *  author: pauletti
 *  date: Aug 21, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_uniform_quad_cache.h>
#include <igatools/basis_functions/bspline_element_accessor.h>

template<int dim, int range, int rank = 1>
void bspline_iterator(const int deg = 1, const int n_knots = 3)
{

    OUTSTART

    auto grid = CartesianGrid<dim>::create(n_knots);
    using Space = BSplineSpace<dim, range, rank>;
    using SpaceCache = BSplineUniformQuadCache<dim, range, rank>;
    auto space = Space::create(deg, grid);

    const int n_qp = 1;
    QGauss< dim > quad(n_qp);

    {
        SpaceCache cache(space, ValueFlags::value, quad);
        auto elem = space->begin();
        cache.init_element_cache(elem);
        for (; elem != space->end(); ++elem)
        {
            cache.fill_element_cache(elem);
            elem->get_basis_values().print_info(out);
        }
    }


    {
        SpaceCache cache(space, ValueFlags::gradient, quad);
        auto elem = space->begin();
        cache.init_element_cache(elem);
        for (; elem != space->end(); ++elem)
        {
            cache.fill_element_cache(elem);
            elem->get_basis_gradients().print_info(out);
        }
    }



    {
        SpaceCache cache(space, ValueFlags::hessian, quad);
        auto elem = space->begin();
        cache.init_element_cache(elem);
        for (; elem != space->end(); ++elem)
        {
            cache.fill_element_cache(elem);
            elem->get_basis_hessians().print_info(out);
        }
    }
}

int main()
{
    out.depth_console(10);

    bspline_iterator<1, 1>();
    bspline_iterator<2, 1>();

    return 0;
}
