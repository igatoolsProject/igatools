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
 *  Test for BSplineBasis using a non range-homogeneous space
 *  Evaluates values gradients and derivatives at two quad point
 *  on each element without the use of the cache.
 *  This test computes the same quantities of bspline_element_iterator_04.cpp
 *
 *  author: martinelli
 *  date: May 08, 2014
 *
 */
// TODO (pauletti, Jun 2, 2014): no need to print info of the spaces
// TODO (pauletti, Jun 2, 2014): rename this file to have 04 at the end
#include "../tests.h"


#include <igatools/basis_functions/bspline_space.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_element_accessor.h>

template <int dim>
shared_ptr<BSplineSpace<dim,dim,1> >
create_space(const int num_knots) ;

template <>
shared_ptr<BSplineSpace<2,2,1> >
create_space<2>(const int num_knots)
{
    auto knots = CartesianGrid<2>::create(num_knots);

    typename BSplineSpace<2,2,1>::DegreeTable degree = { {{3,2}},
        {{2,3}}
    } ;

    return BSplineSpace<2,2,1>::create(degree, knots) ;
}

template <>
shared_ptr<BSplineSpace<3,3,1> >
create_space<3>(const int num_knots)
{
    auto knots = CartesianGrid<3>::create(num_knots);

    typename BSplineSpace<3,3,1>::DegreeTable  degree = { {{3,2,2}},
        {{2,3,2}},
        {{2,2,3}}
    } ;

    return BSplineSpace<3,3,1>::create(degree, knots) ;
}


template< int dim_domain>
void do_test()
{
    out << "domain, range and rank: " << dim_domain << "," << dim_domain << ",1" << endl ;

    const int num_knots = 3 ;

    auto space = create_space<dim_domain>(num_knots) ;

    const int n_points = 2;
    QGauss< dim_domain > quad(n_points) ;
    const auto eval_points = quad.get_points().get_flat_cartesian_product();

    auto elem = space->begin();
    for (; elem != space->end(); ++elem)
    {
        const auto values = elem->evaluate_basis_values_at_points(eval_points);
        out << "Values:" << endl ;
        values.print_info(out);
    }

    {
        auto elem = space->begin();
        for (; elem != space->end(); ++elem)
        {
            const auto gradients = elem->evaluate_basis_gradients_at_points(eval_points);
            out << "Gradients:" << endl ;
            gradients.print_info(out);
        }
    }
    //*/

    {
        auto elem = space->begin();
        for (; elem != space->end(); ++elem)
        {
            const auto hessians = elem->evaluate_basis_hessians_at_points(eval_points);
            out << "Hessians:" << endl ;
            hessians.print_info(out);
        }
    }
    //*/
}


int main()
{
    out.depth_console(10); //to be removed after test finished

    do_test<2>() ;
    do_test<3>() ;

    return 0;
}
