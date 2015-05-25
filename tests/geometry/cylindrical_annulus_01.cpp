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
 *  Test cylindrical annulus mapping values.
 *  author: antolin
 *  date: 2014-03-18
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_lib.h>
#include <igatools/functions/identity_function.h>
#include <igatools/functions/function_element.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>


template <int dim>
void mapping_values()
{
    using Function = functions::CylindricalAnnulus<dim>;

    auto flag = ValueFlags::point | ValueFlags::value |
                ValueFlags::gradient |
                ValueFlags::hessian |
                ValueFlags::measure|
                ValueFlags::w_measure;

    auto quad = QGauss<dim>(1);
    auto grid = CartesianGrid<dim>::create();

    auto F = Function::create(grid, IdentityFunction<dim>::create(grid),
                              1, 2, 0, 2.0, 0.0, numbers::PI/2.0);


    using Mapping   = Mapping<dim, 0>;
    using ElementIt = typename Mapping::ElementIterator;
    auto map = Mapping::create(F);
    map->reset(flag, quad);

//    ElementIt elem(map, 0);
//    ElementIt end(map, IteratorState::pass_the_end);
    ElementIt elem = map->begin();
    ElementIt end = map->end();

    map->template init_cache<dim>(elem);
    for (; elem != end; ++elem)
    {
        map->template fill_cache<dim>(elem, 0);

        out << "Points:" << endl;
        elem->get_points().print_info(out);
        out << endl;
        out << "Values:" << endl;
        elem->template get_values<_Value, dim>(0).print_info(out);
        out << endl;
        out << "Gradients:" << endl;
        elem->template get_values<_Gradient, dim>(0).print_info(out);
        out << endl;
        out << "Hessians:" << endl;
        elem->template get_values<_Hessian, dim>(0).print_info(out);
        out << endl;
        out << "Measure:" << endl;
        elem->template get_measures<dim>(0).print_info(out);
        out << endl;
        out << "weight * measure:" << endl;
        elem->template get_w_measures<dim>(0).print_info(out);
        out << endl;
    }
}


int main()
{
    out.depth_console(10);

    mapping_values<3>();

    return 0;
}


#if 0
int main()
{
    out.depth_console(20);


    auto map = CylindricalAnnulus::create(1, 2, 0, 2.0, 0.0, numbers::PI/2.0);

    const int num_pts = 1 ;
    QGauss<3> quad(num_pts) ;

    auto elem = map->begin();
    ValueFlags flag = ValueFlags::point|ValueFlags::map_gradient|ValueFlags::map_hessian;
    elem->init_cache(flag, quad);
    elem->fill_cache();

    auto values = elem->get_map_values();
    auto gradients = elem->get_map_gradients();
    auto hessians = elem->get_map_hessians();

    out << "Points:" << endl;
    quad.get_points().get_flat_cartesian_product().print_info(out);
    out << endl;

    out << "Values (x1,x2,...):" << endl;
    values.print_info(out);
    out << endl;

    out << "Gradients:" << endl;
    gradients.print_info(out);
    out << endl;

    out << "Hessians:" << endl;
    hessians.print_info(out);
    out << endl;
    out << endl;

//    string filename = "cylindrical_map-" + to_string(3) + "d";
//    Writer<3> writer(map, 4);
//    writer.save(filename);


//TODO(pauletti, Apr 27, 2014): the code below do NOT match the expected test
#if 0
    //----------------------------------------------------------------------------------------------
    typedef PushForward<Transformation::h_grad,3,0> push_fwd_t ;
    auto push_forward = push_fwd_t::create(map);

    PushForwardElementAccessor<push_fwd_t> push_fwd_accessor(push_forward, 0) ;

    push_fwd_accessor.init_cache(ValueFlags::tran_value | ValueFlags::tran_gradient, quad) ;

    push_fwd_accessor.fill_cache() ;

    typedef Values<1,1,1> Value ;
    typedef Derivatives<3,1,1,1> Grad ;

    ValueVector< Value > dummy ;

    ValueVector< Grad > grad_phi0_hat(num_pts) ;
    ValueVector< Grad > grad_phi1_hat(num_pts) ;
    ValueVector< Grad > grad_phi2_hat(num_pts) ;

    for (int iPt = 0 ; iPt < num_pts ; iPt++)
    {
        grad_phi0_hat[iPt][0] = 1.0 ;
        out << grad_phi0_hat[iPt] << endl ;


        grad_phi1_hat[iPt][1] = 1.0 ;
        out << grad_phi1_hat[iPt] << endl ;


        grad_phi2_hat[iPt][2] = 1.0 ;
        out << grad_phi2_hat[iPt] << endl ;
    }

    ValueVector< Grad > grad_phi0(num_pts);
    ValueVector< Grad > grad_phi1(num_pts);
    ValueVector< Grad > grad_phi2(num_pts);


    push_fwd_accessor.transform_gradients<1,1>(dummy,grad_phi0_hat,grad_phi0) ;
    push_fwd_accessor.transform_gradients<1,1>(dummy,grad_phi1_hat,grad_phi1) ;
    push_fwd_accessor.transform_gradients<1,1>(dummy,grad_phi2_hat,grad_phi2) ;

    for (int iPt = 0 ; iPt < num_pts ; iPt++)
    {
        out << grad_phi0[iPt] << endl ;

        out << grad_phi1[iPt] << endl ;

        out << grad_phi2[iPt] << endl ;
    }

#endif

    return (0) ;
}
#endif
