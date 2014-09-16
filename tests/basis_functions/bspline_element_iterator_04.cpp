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
 *  Test for BSplineBasis using a non range-homogeneous space
 *  Evaluates values gradients and derivatives at two quad point
 *  on each element
 *
 *  author: martinelli
 *  date: Oct 04, 2013
 *
 */

#include "../tests.h"


#include <igatools/basis_functions/bspline_space.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_element_accessor.h>

using std::shared_ptr;

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

    typename BSplineSpace<3,3,1>::DegreeTable degree = { {{3,2,2}},
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

    space->print_info(out) ;

    const int n_points = 2;
    QGauss<dim_domain> quad(n_points) ;

    using BSplineCache = BSplineUniformQuadCache<dim_domain,dim_domain>;

    {
        BSplineCache cache(space,ValueFlags::value,quad);

        auto elem = space->begin();
        cache.init_element_cache(elem);

        for (; elem != space->end(); ++elem)
        {
            cache.fill_element_cache(elem);
            out << "Values:" << endl ;
            elem->get_basis_values().print_info(out);
        }
    }

    {
        BSplineCache cache(space,ValueFlags::gradient,quad);

        auto elem = space->begin();
        cache.init_element_cache(elem);

        for (; elem != space->end(); ++elem)
        {
            cache.fill_element_cache(elem);
            out << "Gradients:" << endl ;
            elem->get_basis_gradients().print_info(out);
        }
    }

    {
        BSplineCache cache(space,ValueFlags::hessian,quad);

        auto elem = space->begin();
        cache.init_element_cache(elem);

        for (; elem != space->end(); ++elem)
        {
            cache.fill_element_cache(elem);
            out << "Hessians:" << endl ;
            elem->get_basis_hessians().print_info(out);
        }
    }
}


int main()
{
    out.depth_console(10); //to be removed after test finished

    do_test<2>() ;
    do_test<3>() ;

    return 0;
}
