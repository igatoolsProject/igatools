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
#include <igatools/linear_algebra/epetra_vector.h>
#include <igatools/linear_algebra/dof_tools.h>

using namespace EpetraTools;

template<int dim, int range, int rank =  1>
void evaluate_field(const int deg = 1)
{
  OUTSTART

  auto grid = CartesianGrid<dim>::create();

    using Space = BSplineSpace<dim,range,rank>;

    auto space = Space::create_nonconst(deg, grid);
    const auto  num = space->get_num_basis();
    Epetra_SerialComm comm;
    Epetra_Map map(num, num, 0, comm);

  Vector u(map);
  {
    int id = 0 ;
    u[id++] = 0.0 ;
    u[id++] = 1.0 ;

    u[id++] = 0.0 ;
    u[id++] = 1.0 ;

    u[id++] = 0.0 ;
    u[id++] = 0.0 ;

    u[id++] = 1.0 ;
    u[id++] = 1.0 ;
  }

    auto quad = QGauss<dim>::create(2);
    auto flag = space_element::Flags::value |
                space_element::Flags::gradient |
                space_element::Flags::hessian |
                space_element::Flags::divergence;


    auto elem = space->begin();

    auto elem_handler = space->get_elem_handler();

    elem_handler->template set_flags<dim>(flag);
    elem_handler->init_element_cache(elem,quad);
    elem_handler->fill_element_cache(elem);

    using Elem = typename Space::ElementAccessor;
    using _Value = typename Elem::_Value;
    using _Gradient = typename Elem::_Gradient;
    using _Hessian = typename Elem::_Hessian;
    using _Divergence = typename Elem::_Divergence;

  const auto elem_dofs = elem->get_local_to_global(DofProperties::active);

    const auto &loc_coef = u.get_local_coeffs(elem_dofs);
    out.begin_item("Linear combination of basis values:");
    elem->template linear_combination<_Value,dim>(loc_coef,0,DofProperties::active).print_info(out);
    out.end_item();

    out.begin_item("Linear combination of basis gradients:");
    elem->template linear_combination<_Gradient,dim>(loc_coef,0,DofProperties::active).print_info(out);
    out.end_item();

    out.begin_item("Linear combination of basis hessians:");
    elem->template linear_combination<_Hessian,dim>(loc_coef,0,DofProperties::active).print_info(out);
    out.end_item();

    out.begin_item("Linear combination of basis divergence:");
    elem->template linear_combination<_Divergence,dim>(loc_coef,0,DofProperties::active).print_info(out);
    out.end_item();

  OUTEND
}


int main()
{
  out.depth_console(10);

  evaluate_field<2,2>();
  evaluate_field<3,3>();

  return 0;
}
