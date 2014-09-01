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
#include <igatools/basis_functions/multi_patch_space-template.h>

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

    using iga::vector;

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

    //---------------------------------------------------------------
    const int n_patches = 4;

    using PhysSpacePtr = shared_ptr<PhysicalSpace_t<dim>>;
    using iga::vector;
    vector<PhysSpacePtr> phys_spaces;

    const int num_knots = 2 ;
    out << "Dim: " << dim << endl ;
    const int p = 1 ;

    for (int patch_id = 0 ; patch_id < n_patches ; ++patch_id)
    {
        auto grid = CartesianGrid<dim>::create(num_knots);

        auto ref_space = RefSpace_t<dim>::create(p,grid);

        auto map = create_mapping<dim>(ref_space);

        auto push_fwd = PushForward_t<dim>::create(map);

        phys_spaces.push_back(PhysicalSpace_t<dim>::create(ref_space,push_fwd,patch_id));
    }
    //---------------------------------------------------------------


    //---------------------------------------------------------------
    // adding the patches and the interfaces to the multi-patch space structure
    MultiPatchSpace<PhysicalSpace_t<dim>> multi_patch_space;
    multi_patch_space.patch_insertion_open();

    for (auto phys_space : phys_spaces)
        multi_patch_space.add_patch(phys_space);

    multi_patch_space.patch_insertion_close();

    multi_patch_space.interface_insertion_open();
    multi_patch_space.add_interface(InterfaceType::C0_strong,phys_spaces[0],1,phys_spaces[1],0);
    multi_patch_space.add_interface(InterfaceType::C0_strong,phys_spaces[2],1,phys_spaces[3],0);
    multi_patch_space.add_interface(InterfaceType::C0_strong,phys_spaces[0],2,phys_spaces[2],3);
    multi_patch_space.add_interface(InterfaceType::C0_strong,phys_spaces[1],2,phys_spaces[3],3);
    multi_patch_space.interface_insertion_close();
    //---------------------------------------------------------------

#ifdef USE_GRAPH
    multi_patch_space.build_graph();
#endif

    const auto space_manager = multi_patch_space.get_space_manager();

    //---------------------------------------------------------------
    // adding (manually) the equality constraints that ensures C0 strong continuity
    space_manager->equality_constraints_open();

    // patch 0, side 1 --- patch 1, side 0
    space_manager->add_equality_constraint(1, 8);
    space_manager->add_equality_constraint(3,10);
    space_manager->add_equality_constraint(5,12);
    space_manager->add_equality_constraint(7,14);

    // patch 2, side 1 --- patch 3, side 0
    space_manager->add_equality_constraint(17,24);
    space_manager->add_equality_constraint(19,26);
    space_manager->add_equality_constraint(21,28);
    space_manager->add_equality_constraint(23,30);

    // patch 0, side 3 --- patch 2, side 2
    space_manager->add_equality_constraint(2,16);
    space_manager->add_equality_constraint(3,17);
    space_manager->add_equality_constraint(6,20);
    space_manager->add_equality_constraint(7,21);

    // patch 1, side 3 --- patch 3, side 2
    space_manager->add_equality_constraint(10,24);
    space_manager->add_equality_constraint(11,25);
    space_manager->add_equality_constraint(14,28);
    space_manager->add_equality_constraint(15,29);

    space_manager->remove_equality_constraints_redundancies();

    space_manager->equality_constraints_close();
    //---------------------------------------------------------------




    //---------------------------------------------------------------

//    dofs_manager->print_info(out);
    multi_patch_space.print_info(out);

    for (int patch_id = 0 ; patch_id < n_patches ; ++patch_id)
    {
        out << "The local_dof=3 for the space_id="<< patch_id
            << " corresponds to the global_dof="<< space_manager->get_global_dof(patch_id,3) << endl;
        out << "The local_dof=7 for the space_id="<< patch_id
            << " corresponds to the global_dof="<< space_manager->get_global_dof(patch_id,7) << endl;
    }
    //---------------------------------------------------------------



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
