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
 *  @brief  Test for IgGridFunction class on BSpline of degree 2 on a different grid
 *  and parametrization than the one used in ig_grid_function_bspline_02.cpp
 *
 *  author: martinelli
 *  date: 20 Oct 2015
 *
 */


#include "../tests.h"

#include <igatools/functions/ig_grid_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/functions/function_element.h>



template <int dim, int codim=0>
void ig_grid_function_bspline(const int deg = 2)
{
  const int sub_dim = dim;
  using Basis = BSpline<dim, dim+codim>;
  using Function = IgGridFunction<dim,dim+codim>;



  //----------------------------------------------------------------------------------------------
  int n_knots = 2;
  CartesianProductArray<Real , dim> coord ;
  for (int i = 0; i < dim; ++i)
  {
    SafeSTLVector<Real> tmp_coord;
    for (int j = 0; j < n_knots; ++j)
      tmp_coord.push_back(j);
    coord.copy_data_direction(i,tmp_coord);
  }

  auto grid = Grid<dim>::create(coord);
  auto space = Basis::create(SplineSpace<dim,dim+codim>::create(deg,grid));

  IgCoefficients control_pts;

  if (dim == 2)
  {
    int id = 0 ;

    // x coords
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 1.5 ;
    control_pts[id++] = 1.5 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 2.0 ;
    control_pts[id++] = 2.0 ;


    // y coords
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 1.5 ;
    control_pts[id++] = 1.5 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 2.0 ;
    control_pts[id++] = 2.0 ;
    control_pts[id++] = 0.0 ;
  }
  else if (dim == 3)
  {
    int id = 0 ;

    // x coords
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 1.5 ;
    control_pts[id++] = 1.5 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 2.0 ;
    control_pts[id++] = 2.0 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 1.5 ;
    control_pts[id++] = 1.5 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 2.0 ;
    control_pts[id++] = 2.0 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 1.5 ;
    control_pts[id++] = 1.5 ;

    control_pts[id++] = 0.0 ;
    control_pts[id++] = 2.0 ;
    control_pts[id++] = 2.0 ;


    // y coords
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 1.5 ;
    control_pts[id++] = 1.5 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 2.0 ;
    control_pts[id++] = 2.0 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 1.5 ;
    control_pts[id++] = 1.5 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 2.0 ;
    control_pts[id++] = 2.0 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 1.5 ;
    control_pts[id++] = 1.5 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 2.0 ;
    control_pts[id++] = 2.0 ;
    control_pts[id++] = 0.0 ;

    // z coords
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 0.0 ;
    control_pts[id++] = 0.0 ;

    control_pts[id++] = 0.5 ;
    control_pts[id++] = 0.5 ;
    control_pts[id++] = 0.5 ;
    control_pts[id++] = 0.5 ;
    control_pts[id++] = 0.5 ;
    control_pts[id++] = 0.5 ;
    control_pts[id++] = 0.5 ;
    control_pts[id++] = 0.5 ;
    control_pts[id++] = 0.5 ;

    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;
    control_pts[id++] = 1.0 ;

  }
  auto F = Function::create(space, control_pts);


  auto quad = QGauss<sub_dim>::create(3);
  using Flags = grid_function_element::Flags;
  auto flag =  Flags::D0|
               Flags::D1|
               Flags::D2;

  auto cache_handler = F->create_cache_handler();
  cache_handler->template set_flags<sub_dim>(flag);

  auto elem = F->cbegin();
  auto end  = F->cend();
  const int s_id = 0;

  cache_handler->init_cache(*elem,quad);
  using D0 = grid_function_element::template _D<0>;
  using D1 = grid_function_element::template _D<1>;
  using D2 = grid_function_element::template _D<2>;

  out.begin_item("ig_grid_function_bspline<" + std::to_string(dim) + ">");

  int elem_id = 0;
  for (; elem != end; ++elem, ++elem_id)
  {
    out.begin_item("Element: " + std::to_string(elem_id));

    cache_handler->template fill_cache<sub_dim>(elem, s_id);

    out.begin_item("Values:");
    elem->template get_values_from_cache<D0,sub_dim>(s_id).print_info(out);
    out.end_item();

    out.begin_item("Gradients:");
    elem->template get_values_from_cache<D1,sub_dim>(s_id).print_info(out);
    out.end_item();

    out.begin_item("Hessians:");
    elem->template get_values_from_cache<D2,sub_dim>(s_id).print_info(out);
    out.end_item();

    out.end_item();
  }
  out.end_item();


}

int main()
{
  out.depth_console(10);

  ig_grid_function_bspline<2>();
  ig_grid_function_bspline<3>();

}
