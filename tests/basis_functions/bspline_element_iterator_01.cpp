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
 *  Test for the Bspline space element iterator using
 *  the uniform quad global cache, passing flags one at a time
 *
 *  author: pauletti
 *  date: Aug 28, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/new_bspline_space.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/bspline_element.h>

template<int dim, int range, int rank = 1>
void bspline_iterator(const int deg = 2)
{
    OUTSTART

    auto grid = CartesianGrid<dim>::create();
    using Space = NewBSplineSpace<dim, range, rank>;
    using ElementHandler = typename Space::ElementHandler;
    auto space = Space::create(deg, grid);

    const int n_qp = 3;
    QGauss< dim > quad(n_qp);
    auto flag = NewValueFlags::value|NewValueFlags::gradient|NewValueFlags::hessian;
    ElementHandler cache(space, flag, quad);

    auto elem = space->begin();
    cache.init_element_cache(elem);
    cache.fill_element_cache(elem);

    auto values    = elem->template get_basis_ders<0,0>(0);//
//    auto gradients = elem->template get_basis_ders<0,1>(0);//elem->get_basis_gradients();
//    auto hessians  = elem->template get_basis_ders<0,2>(0);//get_basis_hessians();
//
//    out.begin_item("Values basis functions:");
//    values.print_info(out);
//    out.end_item();
//
//    out.begin_item("Gradients basis functions:");
//    gradients.print_info(out);
//    out.end_item();
//
//    out.begin_item("Hessians basis functions:");
//    hessians.print_info(out);
//    out.end_item();
}


int main()
{
    out.depth_console(10);

    bspline_iterator<1,1>();
//    bspline_iterator<1,2>();
//    bspline_iterator<1,3>();
//
//    bspline_iterator<2,1>();
//    bspline_iterator<2,2>();
//    bspline_iterator<2,3>();
//
//    bspline_iterator<3,1>();
//    bspline_iterator<3,3>();

    return 0;
}
