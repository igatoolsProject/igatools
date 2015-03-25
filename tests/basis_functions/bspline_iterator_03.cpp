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
 *  Test for the iterator evaluate field values and derivatives
 *  author: pauletti
 *  date: 2014-10-23
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/dof_tools.h>

template<LAPack la_pack, int dim, int range, int rank =  1>
void evaluate_field(const int deg = 1)
{
    OUTSTART

    auto grid = CartesianGrid<dim>::create();

    using Space = BSplineSpace<dim,range,rank>;
    using ElementHandler = typename Space::ElementHandler;

    auto space = Space::create(deg, grid);

    Vector<la_pack> u(space->get_num_basis());
    {
        int id = 0 ;
        u(id++) = 0.0 ;
        u(id++) = 1.0 ;

        u(id++) = 0.0 ;
        u(id++) = 1.0 ;

        u(id++) = 0.0 ;
        u(id++) = 0.0 ;

        u(id++) = 1.0 ;
        u(id++) = 1.0 ;
    }

    QGauss<dim> quad(2) ;
    const auto flag = ValueFlags::value|ValueFlags::gradient;
    auto cache1 = ElementHandler::create(space);
    cache1->reset(flag, quad);

    auto elem = space->begin();
    cache1->init_element_cache(elem);
    cache1->fill_element_cache(elem);

    const auto elem_dofs = elem->get_local_to_global(DofProperties::active);
    const auto &loc_coef = u.get_local_coefs(std::set<Index>(elem_dofs.begin(),elem_dofs.end()));
    elem->template linear_combination<0,dim>(loc_coef,0,DofProperties::active).print_info(out);
    out << endl;
    elem->template linear_combination<1,dim>(loc_coef,0,DofProperties::active).print_info(out);
    out << endl;

    OUTEND
}


int main()
{
    out.depth_console(10);

#if defined(USE_TRILINOS)
    const auto la_pack_tpetra = LAPack::trilinos_tpetra;
    out.begin_item("Using Trilinos/TPetra:");
    evaluate_field<la_pack_tpetra,2,2>();
    evaluate_field<la_pack_tpetra,3,3>();
    out.end_item();


    const auto la_pack_epetra = LAPack::trilinos_epetra;
    out.begin_item("Using Trilinos/EPetra:");
    evaluate_field<la_pack_epetra,2,2>();
    evaluate_field<la_pack_epetra,3,3>();
    out.end_item();

#elif defined(USE_PETSC)
    const auto la_pack = LAPack::petsc;

    const auto la_pack_petsc = LAPack::petsc;
    out.begin_item("Using PETSc:");
    evaluate_field<la_pack_petsc,2,2>();
    evaluate_field<la_pack_petsc,3,3>();
    out.end_item();

#endif


    return 0;
}
