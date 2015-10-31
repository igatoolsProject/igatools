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
 *  @brief Mapping using an IgFunction
 *  @author pauletti
 *  @date 2014-10-23
 */

#include "../tests.h"

#include <igatools/functions/ig_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/functions/function_element.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>


template<int dim, int codim, int sub_dim=dim>
void ig_mapping(const int deg = 1)
{
  OUTSTART

  using Space = BSpline<dim, dim+codim>;
  //  using RefSpace = ReferenceSpace<dim, dim+codim>;
  using Function = IgFunction<dim,0,dim+codim,1>;
  using Mapping   = Domain<dim, codim>;


  auto flag =  ValueFlags::value| ValueFlags::gradient
               | ValueFlags::hessian;
  auto quad = QGauss<dim>(2);
  auto grid = Grid<dim>::create(3);

  auto space = Space::create(deg, grid);

  auto c_p = EpetraTools::create_vector(*space, DofProperties::active,Epetra_SerialComm());
  auto &coeff = *c_p;

  coeff[0] = 1.;
  auto F = Function::create(space, c_p);

  auto map = Mapping::create(F);
  map->reset(flag, quad);

  auto elem = map->begin();
  auto end  = map->end();
  const int s_id = 0;

  map->init_cache(elem, Topology<sub_dim>());
  for (; elem != end; ++elem)
  {
    map->fill_cache(elem, Topology<sub_dim>(), s_id);

    elem->template get_values<_Value,sub_dim>(s_id).print_info(out);
    out << endl;
    elem->template get_values<_Gradient,sub_dim>(s_id).print_info(out);
    out << endl;
    elem->template get_values<_Hessian,sub_dim>(s_id).print_info(out);
    out << endl;
  }

  OUTEND
}


int main()
{
  ig_mapping<2,0>();
  ig_mapping<3,0>();

  return 0;
}

