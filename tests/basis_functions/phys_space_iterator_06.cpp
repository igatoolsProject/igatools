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
 *  Test for the evaluation of physical space basis functions
 *  values and gradients with the identity mapping
 *
 *  author: pauletti
 *  date: 2013-10-02
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid_function_lib.h>

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/physical_space_basis.h>
#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/basis_functions/phys_space_element_handler.h>

template<int dim, int codim=0>
auto
create_function(shared_ptr<const Grid<dim>> grid)
{

  using Function = grid_functions::LinearGridFunction<dim,dim+codim>;
  typename Function::Value    b;
  typename Function::Gradient A;

  for (int j=0; j<dim; j++)
    A[j][j] = 1.;

  return Function::const_create(grid,A, b);
}


template <int dim, int k, int range=1, int rank=1, int codim = 0>
void elem_values(const int n_knots = 2, const int deg=1, const int n_qp = 1)
{
  OUTSTART
  using BspSpace = BSpline<dim, range, rank>;
//    using RefSpace = ReferenceSpaceBasis<dim, range,rank>;
  using Basis = PhysicalSpaceBasis<dim,range,rank,codim>;

  auto grid  = Grid<dim>::const_create(n_knots);

  auto ref_space = BspSpace::const_create(
                     SplineSpace<dim,range,rank>::const_create(deg,grid));
  auto map_func = create_function(grid);

  auto space = Basis::const_create(
                 ref_space,
                 Domain<dim,0>::const_create(map_func), Transformation::h_grad);


  auto quad = QGauss<k>::create(n_qp);

  using space_element::Flags;
  auto flag = Flags::value|
              Flags::gradient|
              Flags::hessian |
              Flags::point;

  auto elem_filler = space->create_cache_handler();
  elem_filler->template set_flags<k>(flag);

  auto elem = space->begin();
  auto end = space->end();
  elem_filler->init_face_cache(elem,quad);


  using space_element::_Value;
  using space_element::_Gradient;
  using space_element::_Hessian;
  int elem_id = 0;
  for (; elem != end; ++elem)
  {
    const auto &grid_elem = elem->get_grid_element();
    if (grid_elem.is_boundary())
    {
      out.begin_item("Element " + std::to_string(elem_id));
      out << "Element index: " << grid_elem.get_index() << endl;
      for (auto &s_id : UnitElement<dim>::template elems_ids<k>())
      {
        if (grid_elem.is_boundary(s_id))
        {
          out.begin_item("Face " + std::to_string(s_id));
          elem_filler->fill_face_cache(elem,s_id);

          out.begin_item("Values: ");
          elem->template get_basis_data<_Value, k>(s_id,DofProperties::active).print_info(out);
          out.end_item();

          out.begin_item("Gradients: ");
          elem->template get_basis_data<_Gradient, k>(s_id,DofProperties::active).print_info(out);
          out.end_item();

          out.begin_item("Hessians: ");
          elem->template get_basis_data<_Hessian, k>(s_id,DofProperties::active).print_info(out);
          out.end_item();

          out.end_item();
        }
      }
      out.end_item();
    }
    ++elem_id;
  }
  OUTEND

}


int main()
{
  out.depth_console(10);

  const int p = 1;

  for (int num_knots = 2; num_knots<4; ++num_knots)
  {
    elem_values<2,1>(num_knots, p, 2);
    elem_values<3,2>(num_knots, p, 2);
  }

  return 0;
}




