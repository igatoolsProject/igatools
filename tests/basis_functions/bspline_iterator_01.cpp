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
 *  Test for the Bspline space element iterator using
 *  the uniform quad global cache, passing flags one at a time
 *
 *  author: pauletti
 *  date: Aug 28, 2014
 *
 */
// TODO (pauletti, Oct 23, 2014): this test is very similar bspline_iterator_14
#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/bspline_element.h>

template<int dim, int k=dim, int range = 1, int rank = 1>
void bspline_iterator(const int deg = 2,const int n_qp = 3)
{
  OUTSTART

  auto grid = Grid<dim>::create();
  auto space = SplineSpace<dim,range,rank>::create(deg,grid);
  using Basis = BSpline<dim, range, rank>;
  auto basis = Basis::create(space);


  auto quad = QGauss<k>::create(n_qp);
  auto flag = space_element::Flags::value |
              space_element::Flags::gradient |
              space_element::Flags::hessian;

  auto elem = basis->begin();

  auto elem_handler = basis->create_cache_handler();

  elem_handler->template set_flags<k>(flag);
  elem_handler->template init_cache<k>(*elem,quad);

  using Elem = typename Basis::ElementAccessor;
  using _Value = typename Elem::_Value;
  using _Gradient = typename Elem::_Gradient;
  using _Hessian = typename Elem::_Hessian;

  for (auto &s_id : UnitElement<dim>::template elems_ids<k>())
  {
    elem_handler->template fill_cache<k>(*elem,0);

    out << "Sub Element: " << s_id << endl;
    auto values    = elem->template get_basis<_Value,k>(s_id,DofProperties::active);
    auto gradients = elem->template get_basis<_Gradient,k>(s_id,DofProperties::active);
    auto hessians  = elem->template get_basis<_Hessian,k>(s_id,DofProperties::active);

    out.begin_item("Values basis functions:");
    values.print_info(out);
    out.end_item();

    out.begin_item("Gradients basis functions:");
    gradients.print_info(out);
    out.end_item();

    out.begin_item("Hessians basis functions:");
    hessians.print_info(out);
    out.end_item();
  }
}


template<int dim, int k=dim, int range = 1, int rank = 1>
void bspline_iterator_active_dofs(const int deg = 2,const int n_qp = 3)
{
  OUTSTART

  auto grid = Grid<dim>::create();
  auto space = SplineSpace<dim,range,rank>::create(deg,grid);
  using Basis = BSpline<dim, range, rank>;
  auto basis = Basis::create(space);

  auto dof_distribution = basis->get_ptr_dof_distribution();
  //dof_distribution->add_dofs_property(DofProperties::active);
  for (const auto dof: dof_distribution->get_dofs_view())
    if (dof % 2 == 0)
      dof_distribution->set_dof_property_status(DofProperties::active,dof,true);
    else
      dof_distribution->set_dof_property_status(DofProperties::active,dof,false);

  auto quad = QGauss<k>::create(n_qp);
  auto flag = space_element::Flags::value |
              space_element::Flags::gradient |
              space_element::Flags::hessian;

  auto elem = basis->begin();

  auto elem_handler = basis->create_cache_handler();

  elem_handler->template set_flags<k>(flag);
  elem_handler->init_element_cache(elem,quad);

  using Elem = typename Basis::ElementAccessor;
  using _Value = typename Elem::_Value;
  using _Gradient = typename Elem::_Gradient;
  using _Hessian = typename Elem::_Hessian;

  for (auto &s_id : UnitElement<dim>::template elems_ids<k>())
  {
    elem_handler->fill_element_cache(elem);

    out << "Sub Element: " << s_id << endl;
    auto values    = elem->template get_basis<_Value,k>(s_id,DofProperties::active);
    auto gradients = elem->template get_basis<_Gradient,k>(s_id,DofProperties::active);
    auto hessians  = elem->template get_basis<_Hessian,k>(s_id,DofProperties::active);

    out.begin_item("Values basis functions:");
    values.print_info(out);
    out.end_item();

    out.begin_item("Gradients basis functions:");
    gradients.print_info(out);
    out.end_item();

    out.begin_item("Hessians basis functions:");
    hessians.print_info(out);
    out.end_item();
  }
}


int main()
{
  out.depth_console(0);

  bspline_iterator<1,1,1>();
  bspline_iterator_active_dofs<1,1,1>();

  bspline_iterator<1,1,2>();
  bspline_iterator_active_dofs<1,1,2>();

  bspline_iterator<1,1,3>();
  bspline_iterator_active_dofs<1,1,3>();

  bspline_iterator<2,2,1>();
  bspline_iterator_active_dofs<2,2,1>();

  bspline_iterator<2,2,2>();
  bspline_iterator_active_dofs<2,2,2>();

  bspline_iterator<2,2,3>();
  bspline_iterator_active_dofs<2,2,3>();

  bspline_iterator<3,3,1>();
  bspline_iterator_active_dofs<3,3,1>();

  bspline_iterator<3,3,3>();
  bspline_iterator_active_dofs<3,3,3>();

  return 0;
}
