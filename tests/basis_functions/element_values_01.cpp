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
 *  Test for ElementValues class.
 *  author: pauletti
 *  date: 2013-01-18
 *  updated: 2012-04-02
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/mapping_lib.h>


// output the values of basis functions in the reference domain
template< int dim, int dim_range, int rank >
void do_test(const int p, const int num_knots)
{
    out << "Values, dim: " << dim <<" degree: " << p << endl;

    typedef BSplineSpace< dim, dim_range, rank > space_ref_t;
    typedef PushForward<Transformation::h_grad,dim,0> push_forward_t ;
    typedef PhysicalSpace<space_ref_t,push_forward_t> space_phys_t;

    auto knots = CartesianGrid<dim>::create(num_knots);
    auto space = space_ref_t::create(knots, p);
    auto map = IdentityMapping<dim>::create(knots);
    auto push_forward = push_forward_t::create(map) ;
    auto phys_space = space_phys_t::create(space, push_forward) ;

    QTrapez<dim> quad;
    const int n_qpoints =  quad.get_num_points();

    ValueFlags flags = ValueFlags::value|
                       ValueFlags::gradient |
                       ValueFlags::hessian |
                       ValueFlags::w_measure ;


    const int n_basis = space->get_num_basis_per_element();

    typename space_phys_t::ElementIterator elem = phys_space->begin() ;
    typename space_phys_t::ElementIterator endc = phys_space->end();
    elem->init_values(flags, quad);

    for (; elem != endc; ++elem)
    {
        elem->fill_values();

        out << "Element" << elem->get_flat_index() << endl;

        for (int i=0; i<n_basis; ++i)
        {
            out << "phi_" << i <<"= \n";
            for (int j=0; j<n_qpoints; ++j)
                out << elem->get_basis_value(i,j) << endl;
        }

        for (int i=0; i<n_basis; ++i)
        {
            out << "grad(phi_" << i <<")= \n";
            for (int j=0; j<n_qpoints; ++j)
                out << elem->get_basis_gradient(i,j) << endl;
        }

        for (int i=0; i<n_basis; ++i)
        {
            out << "hessian(phi_" << i <<")= \n";
            for (int j=0; j<n_qpoints; ++j)
                out << elem->get_basis_hessian(i,j) << endl;
        }

        out << "w(qp) * det(DF) = \n";
        for (int j=0; j<n_qpoints; ++j)
            out << elem->get_w_measures()[j] << endl;

        out << endl;
    }
}


//Now a test the values on a deformed domain
template< int dim_ref_domain, int dim_phys_domain, int dim_range, int rank >
void do_test1(const int p)
{
    const int codim = dim_phys_domain - dim_ref_domain;
    out << "Values, dim: " << dim_ref_domain <<" degree: " << p << endl;

    typedef BSplineSpace< dim_ref_domain, dim_range, rank > space_ref_t;
    typedef PushForward<Transformation::h_grad,dim_ref_domain,codim> push_forward_t ;
    typedef PhysicalSpace<space_ref_t,push_forward_t> space_phys_t;


    const int num_knots = 3;
    auto grid = CartesianGrid<dim_ref_domain>::create(num_knots);
    auto space = space_ref_t::create(grid, p);



    Derivatives< dim_ref_domain,dim_phys_domain, 1, 1 > A ;
    Point< dim_phys_domain > b;
    for (int i=0; i<dim_ref_domain; ++i)
        A[i][i] = i+1;
    for (int i=0; i<dim_ref_domain; ++i)
        b[i] = i+1;
    auto map = LinearMapping<dim_ref_domain,codim>::create(grid, A, b);

    auto push_forward = push_forward_t::create(map) ;
    auto phys_space = space_phys_t::create(space, push_forward) ;

    QTrapez<dim_ref_domain> quad;
    const int n_qpoints =  quad.get_num_points();
    ValueFlags flags = ValueFlags::value|
                       ValueFlags::gradient |
                       ValueFlags::hessian |
                       ValueFlags::w_measure ;

    const int n_basis = space->get_num_basis_per_element();

    typename space_phys_t::ElementIterator elem = phys_space->begin() ;
    typename space_phys_t::ElementIterator endc = phys_space->end();
    elem->init_values(flags, quad);
    for (; elem != endc; ++elem)
    {
        elem->fill_values();

        out << "Element: " << elem->get_flat_index() << endl;

        for (int i=0; i<n_basis; ++i)
        {
            out << "phi_" << i <<"= \n";
            for (int j=0; j<n_qpoints; ++j)
                out << elem->get_basis_value(i,j) << endl;
        }

        for (int i=0; i<n_basis; ++i)
        {
            out << "grad(phi_" << i <<")= \n";
            for (int j=0; j<n_qpoints; ++j)
                out << elem->get_basis_gradient(i,j) << endl;
        }

        for (int i=0; i<n_basis; ++i)
        {
            out << "hessian(phi_" << i <<")= \n";
            for (int j=0; j<n_qpoints; ++j)
                out << elem->get_basis_hessian(i,j) << endl;
        }

        out << "w(qp) * det(DF) = \n";
        for (int j=0; j<n_qpoints; ++j)
            out << elem->get_w_measures()[j] << endl;

        out << endl;
    }
}


int main()
{
    out.depth_console(20);

    for (int num_knots = 2; num_knots<4; ++num_knots)
        for (int p=0; p<2; p++)
        {
            do_test<1,1,1>(p, num_knots);
            do_test<2,1,1>(p, num_knots);
            do_test<3,1,1>(p, num_knots);
        }



    for (int p=0; p<2; p++)
    {
        do_test1< 1, 1 ,1, 1>(p) ;
        do_test1< 2, 2, 1, 1>(p) ;
        //  do_test1< 3, 3, 1, 1>(p) ;
    }


    return 0;
}

