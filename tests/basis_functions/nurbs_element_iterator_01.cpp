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

// TODO (pauletti, Jun 10, 2014): write appropriate header comment

#include "../tests.h"

#include <igatools/base/exceptions.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element_accessor.h>



template< int dim, int range, int rank = 1>
void do_test()
{
    const int r = 2;

    out << "do_test<" << dim << "," << range << ">" << endl ;

    using Space = NURBSSpace< dim, range, rank >;
    using WeightsTable = typename Space::WeightsTable;
    using DegreeTable = typename Space::DegreeTable;
    auto  knots = CartesianGrid<dim>::create();

    auto degree = TensorIndex<dim> (r);
    DegreeTable deg(degree);

    auto  bsp = BSplineSpace<dim, range, rank >::create(deg, knots);
    WeightsTable weights;
    const auto n_basis = bsp->get_num_basis_table();
    for (auto comp : Space::components)
        weights(comp).resize(n_basis(comp),1.0);

    auto space = Space::create(deg, knots, weights);

    const int n_points = 3 ;
    QGauss< dim > quad_scheme(n_points) ;

    auto element     = space->begin();
    auto end_element = space->end();


    // initialize the cache of the NURBSSpaceElementAccessor
    element->init_values(ValueFlags::value |
                         ValueFlags::gradient |
                         ValueFlags::hessian,
                         quad_scheme) ;


    out.push("\t") ;
    for (int j = 0 ; element != end_element ; ++j, ++element)
    {
        // fill the cache (with basis functions values, first and second derivatives) at the evaluation points
        element->fill_values() ;


        out << "Element: " << j << endl ;

        out.push("\t") ;


        out << "Values basis functions:" << endl ;
        auto values = element->get_basis_values() ;
        values.print_info(out) ;
        out << endl ;


        out << "Gradients basis functions:" << endl ;
        ValueTable< Derivatives<dim, range, rank, 1 > > gradients = element->get_basis_gradients() ;
        gradients.print_info(out) ;
        out << endl ;


        out << "Hessians basis functions:" << endl ;
        ValueTable< Derivatives<dim, range, rank, 2 > > hessians = element->get_basis_hessians() ;
        hessians.print_info(out) ;
        out << endl ;


        out.pop() ;
    }
    out.pop() ;

    out << endl ;
}


int main()
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
