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

// TODO (pauletti, Jun 11, 2014): put appropriate header or delete

/*
 * Test to figure out gradient bug
 */
#include "../tests.h"

#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/identity_mapping.h>

#include <igatools/linear_algebra/dense_matrix.h>


template<int dim>
void assemble_stiffness_matrix(const int n_knots, const int deg)
{

    auto grid = CartesianGrid<dim>::create(n_knots);
    typedef BSplineSpace<dim> RefSpace;
    auto ref_space = RefSpace::create(deg, grid) ;

    typedef PushForward<Transformation::h_grad,dim,0> PushForward ;
    auto push_forward = PushForward::create(IdentityMapping<dim>::create(grid)) ;

    typedef PhysicalSpace<RefSpace,PushForward> PhysSpace;
    auto phys_space = PhysSpace::create(ref_space, push_forward) ;

    const Quadrature<dim> elem_quad(QGauss<dim>(deg+1)) ;

    ValueFlags flag = ValueFlags::value | ValueFlags::gradient | ValueFlags::w_measure;

    const int n_qpoints =  elem_quad.get_num_points();


    auto elem           = phys_space->begin();
    const auto elem_end = phys_space->end();
    elem->init_cache(flag, elem_quad);

    for (; elem != elem_end; ++elem)
    {
        elem->fill_cache();

        const int n_basis = elem->get_num_basis();

        DenseMatrix loc_mat(n_basis,n_basis);
        loc_mat = 0.0;

        const auto w_meas = elem->get_w_measures();

        for (int i=0; i<n_basis; ++i)
        {
            const auto grad_i = elem->get_basis_gradients(i);

            for (int j=0; j<n_basis; ++j)
            {
                const auto grad_j = elem->get_basis_gradients(j);

                for (int qp = 0; qp < n_qpoints; ++qp)
                    loc_mat(i,j) += scalar_product(grad_i[qp], grad_j[qp]) * w_meas[qp];
            }
        }
        out << loc_mat << endl;
    }

}




int main()
{
    out.depth_console(1);
    const int n_knots = 6;
    const int deg = 1;

    assemble_stiffness_matrix<1>(n_knots, deg);

    return  0;
}
