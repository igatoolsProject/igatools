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
 *  Test for the multi-patch space using some IgMapping.
 *
 *  author: martinelli
 *  date: 2014-05-22
 *
 */

#include "../tests.h"

#include <igatools/geometry/ig_mapping.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/multi_patch_space.h>

template <int dim>
using RefSpace_t = BSplineSpace<dim,dim,1>  ;

template <int dim>
using PushForward_t = PushForward< Transformation::h_grad,dim,0> ;


template <int dim>
using PhysicalSpace_t = PhysicalSpace< RefSpace_t<dim>, PushForward_t<dim> > ;



template <int dim>
shared_ptr< IgMapping< RefSpace_t<dim> > >
create_mapping(shared_ptr<RefSpace_t<dim>> bspline_space)
{
    // bspline_space->print_info(out) ;

    vector<Real> control_pts(bspline_space->get_num_basis());

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

    return shared_ptr< IgMapping< RefSpace_t<dim> > >(new IgMapping< RefSpace_t<dim> >(bspline_space, control_pts)) ;
}



template <int dim>
void test_evaluate()
{
    out << "========== Test dim=" << dim << " --- begin ==========" << endl;


    const int num_knots = 2 ;
    out << "Dim: " << dim << endl ;
    const int p = 1 ;
    auto knots = CartesianGrid<dim>::create(num_knots);

    auto ref_space_0 = RefSpace_t<dim>::create(knots, p);
    auto ref_space_1 = RefSpace_t<dim>::create(knots, p);

    auto map_0 = create_mapping<dim>(ref_space_0);
    auto map_1 = create_mapping<dim>(ref_space_1);
    auto push_forward_0 = PushForward_t<dim>::create(map_0);
    auto push_forward_1 = PushForward_t<dim>::create(map_1);
    // push_forward->print_info(out) ;

//    auto ref_space1 = RefSpace_t<dim>::create(knots, p);
    auto physical_space_0 = PhysicalSpace_t<dim>::create(ref_space_0, push_forward_0, 0) ;
//    physical_space->print_info(out) ;
    auto physical_space_1 = PhysicalSpace_t<dim>::create(ref_space_1, push_forward_1, 1) ;


    MultiPatchSpace<PhysicalSpace_t<dim>> multi_patch_space;
    multi_patch_space.arrangement_open();

    multi_patch_space.add_patch(physical_space_0);
    multi_patch_space.add_patch(physical_space_1);

    multi_patch_space.add_interface(
        InterfaceType::C0_strong,
        physical_space_0,1,
        physical_space_1,0);

    multi_patch_space.arrangement_close();
    multi_patch_space.print_info(out);


    const auto &dofs_manager = multi_patch_space.get_dofs_manager();

    out << "The local_dof=3 for the space_id=0 corresponds to the global_dof="<< dofs_manager.get_global_dof(0,3) << endl;
    out << "The local_dof=7 for the space_id=0 corresponds to the global_dof="<< dofs_manager.get_global_dof(0,7) << endl;
    out << "The local_dof=3 for the space_id=1 corresponds to the global_dof="<< dofs_manager.get_global_dof(1,3) << endl;
    out << "The local_dof=7 for the space_id=1 corresponds to the global_dof="<< dofs_manager.get_global_dof(1,7) << endl;

    out << "========== Test dim=" << dim << " --- end ==========" << endl;
    out << endl << endl;
}

int main()
{
    out.depth_console(10);

//    test_evaluate<1>();
    test_evaluate<2>();
//    test_evaluate<3>();

}
