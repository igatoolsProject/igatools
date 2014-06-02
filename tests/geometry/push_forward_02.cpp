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
 *  Solve a Poisson-like problem.
 *  author: martinelli
 *  date: 2013-01-18
 *
 */
//TODO(pauletti, Apr 27, 2014): header comment is wront
//TODO(pauletti, Apr 27, 2014): seem the same as cylindrical anulus test

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/geometry/push_forward_element_accessor.h>
#include <igatools/io/writer.h>

int main()
{
    out.depth_console(20);


    auto map = CylindricalAnnulus::create(1, 2, 0, 2.0, 0.0, numbers::PI / 2.0);

    const int num_pts = 1 ;
    QGauss<3> quad(num_pts) ;

    auto elem = map->begin();
    ValueFlags flag = ValueFlags::point|ValueFlags::map_gradient|ValueFlags::map_hessian;
    elem->init_values(flag, quad);
    elem->fill_values();

    auto values = elem->get_map_values();
    auto gradients = elem->get_map_gradients();
    auto hessians = elem->get_map_hessians();

    out << "Points:" << endl;
    out << quad.get_points().get_flat_cartesian_product() << endl;
    out << "Values (x1,x2,...):" << endl;
    values.print_info(out);
    gradients.print_info(out);
    hessians.print_info(out);

    string filename = "cylindrical_map-" + to_string(3) + "d";
    Writer<3> writer(map, 4);
    writer.save(filename);




    //----------------------------------------------------------------------------------------------
    typedef PushForward<Transformation::h_grad,3,0> push_fwd_t ;
    auto push_forward = push_fwd_t::create(map);

    PushForwardElementAccessor<push_fwd_t> push_fwd_accessor(push_forward, 0) ;

    push_fwd_accessor.init_values(ValueFlags::tran_value | ValueFlags::tran_gradient, quad) ;

    push_fwd_accessor.fill_values() ;

    typedef Values<1,1> Value_t ;
    typedef Derivatives<3,1,1,1> Grad_t ;

    ValueVector< Value_t > dummy ;

    ValueVector< Grad_t > grad_phi0_hat(num_pts) ;
    ValueVector< Grad_t > grad_phi1_hat(num_pts) ;
    ValueVector< Grad_t > grad_phi2_hat(num_pts) ;

    for (int iPt = 0 ; iPt < num_pts ; iPt++)
    {
        grad_phi0_hat[iPt][0] = 1.0 ;
        out << grad_phi0_hat[iPt] << endl ;


        grad_phi1_hat[iPt][1] = 1.0 ;
        out << grad_phi1_hat[iPt] << endl ;


        grad_phi2_hat[iPt][2] = 1.0 ;
        out << grad_phi2_hat[iPt] << endl ;
    }

    ValueVector< Grad_t > grad_phi0(num_pts);
    ValueVector< Grad_t > grad_phi1(num_pts);
    ValueVector< Grad_t > grad_phi2(num_pts);


    push_fwd_accessor.transform_gradients<1,1>(dummy,grad_phi0_hat,grad_phi0) ;
    push_fwd_accessor.transform_gradients<1,1>(dummy,grad_phi1_hat,grad_phi1) ;
    push_fwd_accessor.transform_gradients<1,1>(dummy,grad_phi2_hat,grad_phi2) ;

    for (int iPt = 0 ; iPt < num_pts ; iPt++)
    {
        out << grad_phi0[iPt] << endl ;

        out << grad_phi1[iPt] << endl ;

        out << grad_phi2[iPt] << endl ;
    }



    return (0) ;
}

