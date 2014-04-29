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
#include <igatools/geometry/identity_mapping.h>

template <int dim>
using RefSpace_t = NURBSSpace<dim>  ;

template <int dim>
using PushForward_t = PushForward<Transformation::h_grad,dim,0> ;

template <int dim>
using PhysicalSpace_t = PhysicalSpace< RefSpace_t<dim>, PushForward_t<dim> > ;


template <class T, int dim>
using ComponentTable = StaticMultiArray<T, RefSpace_t<dim>::range, RefSpace_t<dim>::rank >;


template <int dim>
void test_evaluate()
{
    const int deg = 2;
    auto grid = CartesianGrid<dim>::create();
    grid->refine();
    out << endl;

    auto map = IdentityMapping<dim>::create(grid);
//   map->refine();
//    map->get_grid()->print_info(out);
//    out << endl;



    auto push_forward = PushForward<Transformation::h_grad,dim,0>::create(map);

    auto ref_space = RefSpace_t<dim>::create(grid, deg);

    TensorSize<dim> n_weights_comp;
    for (Index dir_id = 0 ; dir_id < dim ; ++dir_id)
        n_weights_comp(dir_id) = ref_space->get_num_basis(0,dir_id);

    DynamicMultiArray<Real,dim> weights_comp(n_weights_comp);
    if (dim == 1)
    {
        Index id = 0;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;
    }
    else if (dim == 2)
    {
        Index id = 0;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;
    }
    else if (dim == 3)
    {
        Index id = 0;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;

        weights_comp(id++) = 1.0 ;
        weights_comp(id++) = 0.4 ;
        weights_comp(id++) = 0.65 ;
        weights_comp(id++) = 1.0 ;
    }

    ref_space->reset_weights(ComponentTable<DynamicMultiArray<Real,dim>,dim>(weights_comp));

    auto phys_space = PhysicalSpace_t<dim>::create(ref_space, push_forward);


    out << endl;
    out << endl;

    out << "===============================================================" << endl;
    out << "O R I G I N A L     S P A C E" << endl;
    phys_space->print_info(out);
    out << "===============================================================" << endl;
    out << endl;

    out << "===============================================================" << endl;
    out << "R E F I N E D     S P A C E" << endl;
    phys_space->refine_h();
    phys_space->print_info(out);
    out << "===============================================================" << endl;
    out << endl;
}

int main()
{
    out.depth_console(10);

    test_evaluate<1>();
    out << endl ;

    test_evaluate<2>();
    out << endl ;

    test_evaluate<3>();
    out << endl ;

    return 0;
}
