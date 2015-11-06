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

/**
 *  @file
 *  @brief Domain using an IgGridFunction withouth the use of the cache
 *  @author martinelli
 *  @date Nov 06, 2015
 */

#include "../tests.h"

#include <igatools/functions/ig_grid_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/geometry/domain_handler.h>


template<int dim, int codim>
void ig_mapping(const int deg = 1)
{
  OUTSTART
  out.begin_item("ig_mapping<dim=" + std::to_string(dim) +",codim=" + std::to_string(codim) +">");


  auto grid = Grid<dim>::const_create(3);
  auto space = SplineSpace<dim, dim+codim>::const_create(deg, grid);
  auto basis = BSpline<dim, dim+codim>::const_create(space);

  auto c_p = EpetraTools::create_vector(*basis, DofProperties::active,Epetra_SerialComm());
  (*c_p)[0] = 1.;
  auto F = IgGridFunction<dim,dim+codim>::const_create(basis, *c_p);
  auto domain = Domain<dim, codim>::const_create(F);



  auto elem = domain->begin();
  auto end  = domain->end();

  auto quad = QGauss<dim>::create(2);
  int elem_id = 0;
  for (; elem != end; ++elem, ++elem_id)
  {
    out.begin_item("Element " +std::to_string(elem_id));

    out << "Element ID: " << elem->get_index() << std::endl;

    out.begin_item("Points:");
    elem->template evaluate_at_points<domain_element::_Point>(quad).print_info(out);
    out.end_item();

    out.begin_item("Jacobians:");
    elem->template evaluate_at_points<domain_element::_Jacobian>(quad).print_info(out);
    out.end_item();

    out.begin_item("Hessians:");
    elem->template evaluate_at_points<domain_element::_Hessian>(quad).print_info(out);
    out.end_item();

    out.end_item();
  }

  out.end_item();
  OUTEND
}


int main()
{
  ig_mapping<1,0>();
  ig_mapping<2,0>();
  ig_mapping<3,0>();

  ig_mapping<1,1>();
  ig_mapping<2,1>();

  ig_mapping<1,2>();

  return 0;
}
/*
template<int dim>
void ig_mapping(const int deg = 1)
{
  OUTSTART

  using Basis = BSpline<dim,dim>;
  using RefSpace = ReferenceSpaceBasis<dim, dim>;
  using Function = IgFunction<dim,0,dim,1>;
  using Mapping   = Mapping<dim,0>;


  auto quad = QGauss<dim>(2);
  auto grid = Grid<dim>::create(3);

  auto space = Basis::create(deg, grid);
  auto coeff = EpetraTools::create_vector(*space, DofProperties::active,Epetra_SerialComm());
  (*coeff)[0] = 1.;
  auto F = Function::create(space, coeff);

  auto map = Mapping::create(F);

  auto elem = map->begin();
  auto end  = map->end();

  for (; elem != end; ++elem)
  {
    elem->template evaluate_at_points<_Value>(quad).print_info(out);
    out << endl;
    elem->template evaluate_at_points<_Gradient>(quad).print_info(out);
    out << endl;
    elem->template evaluate_at_points<_Hessian>(quad).print_info(out);
    out << endl;
  }

  OUTEND
}


int main()
{
  ig_mapping<2>();
  ig_mapping<3>();

  return 0;
}
//*/
