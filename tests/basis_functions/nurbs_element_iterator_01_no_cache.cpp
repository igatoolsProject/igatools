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

#include "../tests.h"

#include <igatools/base/exceptions.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element_accessor.h>



template< int dim_domain, int dim_range >
void do_test()
{
    out << "do_test<" << dim_domain << "," << dim_range << ">" << endl ;

    auto knots = CartesianGrid<dim_domain>::create();


    // and here we build the NURBSSpace
    const int rank = 1;

    typedef NURBSSpace< dim_domain, dim_range, rank > Space_t ;
    auto space = Space_t::create(knots, 2) ;

    space->print_info(out) ;
    out << endl;
    //----------------------------------------------------------------------------------------------


    //----------------------------------------------------------------------------------------------
    // for the basis functions evaluation we need a set of points (with tensor product structure)
    // to do so, we get the points from a Gauss quadrature scheme with 3 points

    const int n_points = 3 ;
    QGauss< dim_domain > quad_scheme(n_points) ;
    const auto points = quad_scheme.get_points().get_flat_cartesian_product();


    auto element     = space->begin();
    auto end_element = space->end();




    out.push("\t") ;
    for (int j = 0 ; element != end_element ; ++j, ++element)
    {
        out << "Element: " << j << endl ;

        out.push("\t") ;


        out << "Values basis functions:" << endl ;
        const auto values = element->evaluate_basis_values_at_points(points);
        values.print_info(out) ;
        out << endl ;


        out << "Gradients basis functions:" << endl ;
        const auto gradients = element->evaluate_basis_gradients_at_points(points);
        gradients.print_info(out) ;
        out << endl ;


        out << "Hessians basis functions:" << endl ;
        const auto hessians = element->evaluate_basis_hessians_at_points(points);
        hessians.print_info(out) ;
        out << endl ;


        out.pop() ;
    }
    out.pop() ;

    out << endl ;
}


int main(int argc, char *argv[])
{
    do_test< 1, 1 >() ;

    do_test< 1, 2 >() ;

    do_test< 1, 3 >() ;

    do_test< 2, 1 >() ;

    do_test< 2, 2 >() ;
    do_test< 2, 3 >() ;
    do_test< 3, 1 >() ;
//  do_test< 3, 2 >() ;
    do_test< 3, 3 >() ;
//*/
    return (0) ;
}
