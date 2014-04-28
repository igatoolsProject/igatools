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

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/nurbs_element_accessor.h>



template< int dim_domain, int dim_range >
void do_test()
{
    Assert(dim_domain == 1 || dim_domain == 2 || dim_domain == 3, ExcIndexRange(dim_domain, 1, 4)) ;

    out << "do_test<" << dim_domain << "," << dim_range << ">" << endl ;

    auto knots = CartesianGrid<dim_domain>::create(3);

    // and here we build the NURBSSpace
    const int rank = 1;



    typedef NURBSSpace< dim_domain, dim_range, rank > Space_t ;
    auto space = Space_t::create(knots, 2) ;


    // defining and assigning some weights (different from 1.0) to the NURBSSpace
    const auto n_dofs = space->get_num_dofs() ;

    const auto n_dofs_component = n_dofs(0);
    DynamicMultiArray<iga::Real,dim_domain> weights_component(n_dofs_component) ;
    const int n_entries_component = weights_component.flat_size();
    for (int i = 0 ; i < n_entries_component ; ++i)
        weights_component(i) = (i+1) * (1.0 / n_entries_component) ;

    StaticMultiArray<DynamicMultiArray<iga::Real,dim_domain>,dim_range,rank> weights(weights_component) ;

    space->reset_weights(weights) ;

    space->print_info(out) ;
    out << endl;
    //----------------------------------------------------------------------------------------------


    //----------------------------------------------------------------------------------------------
    // for the basis functions evaluation we need a set of points (with tensor product structure)
    // to do so, we get the points from a Gauss quadrature scheme with 3 points

    const int n_points = 3 ;
    QGauss< dim_domain > quad_scheme(n_points) ;
    //----------------------------------------------------------------------------------------------





    //----------------------------------------------------------------------------------------------
    auto element     = space->begin();
    auto end_element = space->end();

    element->init_values(ValueFlags::value|ValueFlags::gradient|ValueFlags::hessian,
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
        auto gradients = element->get_basis_gradients() ;
        gradients.print_info(out) ;
        out << endl ;


        out << "Hessians basis functions:" << endl ;
        auto hessians = element->get_basis_hessians() ;
        hessians.print_info(out) ;
        out << endl ;


        //-------------------------------------------------------------------------------------
        // here we print the sum of all the function values at the different points (they should be always 1.0)
        out << endl ;

        const int n_points = values.get_num_points() ;
        const int n_functions = values.get_num_functions() ;

        vector<iga::Real> sum_values_same_point(n_points,0.0) ;

        for (int ifunc = 0 ; ifunc < n_functions ; ++ifunc)
        {
            const auto values_function = values.get_function_view(ifunc) ;

            for (int jpt = 0 ; jpt < n_points ; ++jpt)
                sum_values_same_point[jpt] += values_function[jpt][0][0] ;
        }

        for (int jpt = 0 ; jpt < n_points ; ++jpt)
            out << "Function values sum at point " << jpt << " = " << sum_values_same_point[jpt] << endl ;

        out << endl ;
        //-------------------------------------------------------------------------------------

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
    do_test< 3, 3 >() ;

    return 0;
}
