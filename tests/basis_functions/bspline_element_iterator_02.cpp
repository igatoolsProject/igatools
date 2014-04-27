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
 *  Test for BSplineBasis
 *  Evaluates values gradients and derivatives at one quad point
 *  on each element
 *
 *  author: pauletti
 *  date: Aug 21, 2013
 *
 */

#include "../tests.h"


#include <igatools/basis_functions/bspline_space.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_element_accessor.h>

template< int dim_domain, int dim_range, int rank >
void do_test(const int degree)
{
    out << "domain, range and rank: " << dim_domain << "," << dim_range << "," << rank << endl ;

    auto knots = CartesianGrid<dim_domain>::create(3);

    typedef BSplineSpace< dim_domain, dim_range, rank > Space_t ;
    auto space = Space_t::create(knots, degree) ;

    const int n_points = 1;
    QGauss< dim_domain > quad(n_points) ;

    auto elem = space->begin();
    elem->init_values(ValueFlags::value, quad);

    for (; elem != space->end(); ++elem)
    {
        elem->fill_values();
        elem->get_basis_values().print_info(out);
    }

    elem = space->begin();
    elem->init_values(ValueFlags::gradient, quad) ;

    for (; elem != space->end(); ++elem)
    {
        elem->fill_values();
        elem->get_basis_gradients().print_info(out);
    }

    elem = space->begin();
    elem->init_values(ValueFlags::hessian, quad) ;

    for (; elem != space->end(); ++elem)
    {
        elem->fill_values();
        elem->get_basis_hessians().print_info(out);
    }

}


int main()
{
    out.depth_console(10); //to be removed after test finished

    do_test< 1, 1, 1 >(1) ;
    do_test< 2, 1, 1 >(1) ;

    return 0;
}
