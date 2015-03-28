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
 *
 *  Test for the evaluation of physical space basis functions
 *  values and gradients with an ig function for the map
 *  author:
 *  date: 2014-11-08
 *
 */

#include "../tests.h"

#include <igatools/base/ig_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/basis_functions/phys_space_element_handler.h>


template<int dim, int codim=0>
auto
create_function(shared_ptr<BSplineSpace<dim, dim + codim>> space)
{
    using Space = ReferenceSpace<dim, dim + codim>;
    using Function = IgFunction<Space>;

    vector<Real> control_pts(space->get_num_basis());
    if (dim == 1)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;
    }
    else if (dim == 2)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
    }
    else if (dim == 3)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

    }

    return Function::create(space, IgCoefficients(*space,DofProperties::active,control_pts));
}



template <int dim, int order = 0, int range=dim, int rank=1, int codim = 0>
void elem_values(const int n_knots = 2, const int deg=1)
{
    OUTSTART
    const int k = dim;
    using BspSpace = BSplineSpace<dim, range, rank>;
    using RefSpace = ReferenceSpace<dim, range,rank>;
    using Space = PhysicalSpace<dim,range,rank,codim, Transformation::h_grad>;
    using ElementHandler = typename Space::ElementHandler;

    auto grid  = CartesianGrid<dim>::create(n_knots);

    auto ref_space = BspSpace::create(deg, grid);
    auto map_func = create_function(ref_space);

    auto space = Space::create(ref_space, map_func);

    const int n_qp = 3;
    auto quad = QGauss<k>(n_qp);

    auto flag = ValueFlags::value |
                ValueFlags::gradient |
                ValueFlags::hessian |
                ValueFlags::divergence |
                ValueFlags::w_measure;

    ElementHandler sp_values(space);
    sp_values.template reset<k> (flag, quad);

    auto elem = space->begin();
    auto end = space->end();
    sp_values.init_element_cache(elem);
    for (; elem != end; ++elem)
    {
        sp_values.fill_element_cache(elem);
        out.begin_item("Element " + std::to_string(elem->get_flat_index()));
        elem->print_info(out);

        out.begin_item("Values: ");
        elem->template get_basis<_Value, k>(0,DofProperties::active).print_info(out);
        out.end_item();

        out.begin_item("Gradients: ");
        elem->template get_basis<_Gradient, k>(0,DofProperties::active).print_info(out);
        out.end_item();

        out.begin_item("Hessians: ");
        elem->template get_basis<_Hessian, k>(0,DofProperties::active).print_info(out);
        out.end_item();

        out.begin_item("Divergences: ");
        elem->template get_basis<_Divergence,k>(0,DofProperties::active).print_info(out);
        out.end_item();

        out.begin_item("W * Measures: ");
        elem->template get_w_measures<k>(0).print_info(out);
        out.end_item();

        out.end_item();
    }
    OUTEND
}

int main()
{
    out.depth_console(10);

    elem_values<1>();
    elem_values<2>();
    elem_values<3>();

}
