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
 *  Test for the physical space element iterator
 *
 *  author: pauletti
 *  date: 2013-10-02
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/geometry/ig_mapping.h>

template <int dim>
using ReferenceSpace = NURBSSpace<dim,dim>;

template <int dim>
using PushFwd = PushForward<Transformation::h_grad,dim,0> ;

template <int dim>
using PhysSpace = PhysicalSpace< ReferenceSpace<dim>, PushFwd<dim> > ;

template <class T, int dim>
using ComponentTable = StaticMultiArray<T, ReferenceSpace<dim>::range, ReferenceSpace<dim>::rank >;

template <int dim>
void test_evaluate()
{

    auto grid = CartesianGrid<dim>::create(2);
    grid->refine();


    const int deg = 2;
    TensorSize<dim> n_weights_dir(deg+2);

    TensorIndex<dim> component_map;
    for (int i = 0 ; i < dim ; ++i)
        component_map[i] = i;

    typename NURBSSpace<dim,dim>::WeightsTable weights(component_map);
    for (auto &weights_comp : weights)
    {
        weights_comp.resize(n_weights_dir);

        Index id = 0;
        for (Index i=0 ; i < pow(4,dim-1) ; ++i)
        {
            weights_comp[id++] = 1.0 ;
            weights_comp[id++] = 0.853553390593274 ;
            weights_comp[id++] = 0.853553390593274 ;
            weights_comp[id++] = 1.0 ;
        }
    }


    auto ref_space = ReferenceSpace<dim>::create(deg,grid,weights);

    vector<Real> control_pts(ref_space->get_num_basis());

    if (dim == 1)
    {
        control_pts[0] = 1.0;
        control_pts[1] = 1.0;
        control_pts[2] = 0.414213562373095;
        control_pts[3] = 0.0;
    }
    else if (dim == 2)
    {
        // 1st comp - 1st row
        control_pts[0] = 1.0;
        control_pts[1] = 1.0;
        control_pts[2] = 0.414213562373095;
        control_pts[3] = 0.0;

        // 1st comp - 2nd row
        control_pts[4] = 1.375;
        control_pts[5] = 1.375;
        control_pts[6] = 0.569543648263006;
        control_pts[7] = 0.0;

        // 1st comp - 3rd row
        control_pts[8] = 2.125;
        control_pts[9] = 2.125;
        control_pts[10] = 0.880203820042827;
        control_pts[11] = 0.0;

        // 1st comp - 4th row
        control_pts[12] = 2.5;
        control_pts[13] = 2.5;
        control_pts[14] = 1.03553390593274;
        control_pts[15] = 0.0;

        // 2nd comp - 1st row
        control_pts[16] = 0.0;
        control_pts[17] = 0.414213562373095;
        control_pts[18] = 1.0;
        control_pts[19] = 1.0;

        // 2nd comp - 2nd row
        control_pts[20] = 0.0;
        control_pts[21] = 0.569543648263006;
        control_pts[22] = 1.375;
        control_pts[23] = 1.375;

        // 2nd comp - 3rd row
        control_pts[24] = 0.0;
        control_pts[25] = 0.880203820042827;
        control_pts[26] = 2.125;
        control_pts[27] = 2.125;

        // 2nd comp - 4th row
        control_pts[28] = 0.0;
        control_pts[29] = 1.035533905932738;
        control_pts[30] = 2.5;
        control_pts[31] = 2.5;
    }
    else if (dim == 3)
    {
        Assert(false,ExcNotImplemented());
        AssertThrow(false,ExcNotImplemented());
    }

    auto map = IgMapping<ReferenceSpace<dim>>::create(ref_space,control_pts);

    auto push_fwd = PushFwd<dim>::create(map);

    auto phys_space = PhysSpace<dim>::create(ref_space,push_fwd);


    out << endl;
    out << endl;

    out << "===============================================================" << endl;
    out << "O R I G I N A L     S P A C E" << endl;
    phys_space->print_info(out);
    out << "===============================================================" << endl;
    out << endl;

    out << "===============================================================" << endl;
    out << "R E F I N E D     S P A C E (2 elements each old element)" << endl;
    phys_space->refine_h();
    phys_space->print_info(out);
    out << "===============================================================" << endl;
    out << endl;


    out << "===============================================================" << endl;
    out << "R E F I N E D     S P A C E (6 elements each old element)" << endl;
    phys_space->refine_h(3);
    phys_space->print_info(out);
    out << "===============================================================" << endl;
    out << endl;

}

int main()
{
    out.depth_console(10);

//    test_evaluate<1>();
    test_evaluate<2>();
//    test_evaluate<3>();

    return 0;
}
