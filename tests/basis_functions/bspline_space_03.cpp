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
 * todo: Header please
 * author: Max?
 */
#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/base/exceptions.h>
#include <igatools/basis_functions/bspline_element_accessor.h>


template< int dim_domain, int dim_range, int rank >
void do_test(const int p)
{

    out << "Domain dim: " << dim_domain;
    out << " Range dim: " << dim_range;
    out << " Degree: " << p <<endl;

    const int n_knots = 4;

    CartesianProductArray< Real, dim_domain > coord(n_knots) ;

    for (int i = 0 ; i < dim_domain ; ++i)
        for (int j = 0 ; j < n_knots ; ++j)
            coord.entry(i,j) = j;

    auto knots = CartesianGrid<dim_domain>::create(coord);
    auto bspline_space = BSplineSpace< dim_domain, dim_range, rank>::create(p, knots) ;

    typename BSplineSpace< dim_domain, dim_range, rank >::ElementIterator
    element = bspline_space->begin(), endc = bspline_space->end();
    for (; element != endc; ++element)
    {
        out << "Element index: " << element->get_flat_index() << endl;
        out << "Global dofs: ";
        element->get_local_to_global().print_info(out);
        out << endl;
    }

    out << endl;

}


int main(int argc, char *argv[])
{
    out.depth_console(20);
    for (int p=0; p<3; p++)
    {
        do_test< 1, 1, 1 >(p) ;
        do_test< 1, 2, 1 >(p) ;
        do_test< 1, 3, 1 >(p) ;
    }

    return (0) ;
}
