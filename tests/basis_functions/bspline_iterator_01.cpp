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
// TODO (pauletti, Oct 23, 2014): this test is very similar bspline_iterator_14
#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/new_bspline_space.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/bspline_element.h>

template<int dim, int k=dim, int range = 1, int rank = 1>
void bspline_iterator(const int deg = 2,const int n_qp = 3)
{
    OUTSTART

    auto grid = CartesianGrid<dim>::create();
    using Space = NewBSplineSpace<dim, range, rank>;
    using ElementHandler = typename Space::ElementHandler;
    auto space = Space::create(deg, grid);


    QGauss<k> quad(n_qp);
    auto flag = NewValueFlags::value|NewValueFlags::gradient
                |NewValueFlags::hessian;
    ElementHandler cache(space);
    cache.template reset<k>(flag, quad);

    auto elem = space->begin();
    cache.template init_cache<k>(elem);

    for (auto &s_id : UnitElement<dim>::template elems_ids<k>())
    {
    	cache.template fill_cache<k>(elem, s_id);

    	out << "Sub Element: " << s_id << endl;
    	auto values    = elem->template get_values<0,k>(s_id);
    	auto gradients = elem->template get_values<1,k>(s_id);
    	auto hessians  = elem->template get_values<2,k>(s_id);

    	out.begin_item("Values basis functions:");
    	values.print_info(out);
    	out.end_item();

    	out.begin_item("Gradients basis functions:");
    	gradients.print_info(out);
    	out.end_item();

    	out.begin_item("Hessians basis functions:");
    	hessians.print_info(out);
    	out.end_item();
    }
}


int main()
{
    out.depth_console(10);

    bspline_iterator<1,1,1>();
    bspline_iterator<1,1,2>();
    bspline_iterator<1,1,3>();

    bspline_iterator<2,2,1>();
    bspline_iterator<2,2,2>();
    bspline_iterator<2,2,3>();

    bspline_iterator<3,3,1>();
    bspline_iterator<3,3,3>();

    return 0;
}
