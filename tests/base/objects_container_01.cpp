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
 *  @brief  Test for Objects container
 *  @author P. Antolin
 *  @date 2015
 */

#include "../tests.h"

#include <igatools/base/objects_container.h>
#include <igatools/geometry/grid.h>

#include <igatools/functions/grid_function_lib.h>
#include <igatools/functions/function_lib.h>
#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/nurbs_element.h>

template< int dim, int range, int rank = 1>
void insert_objects(const std::shared_ptr<ObjectsContainer> container)
{
  using iga::SafeSTLVector;
  SafeSTLVector<Real> coord_x {0,1,2,3,4};
  SafeSTLVector<Real> coord_y {5,6,7,8};
  SafeSTLVector<Real> coord_z {9, 10, 11};

  SafeSTLArray<SafeSTLVector<Real>, dim> coord;
  TensorIndex<dim> degree;

  if (dim == 1)
  {
    coord[0] = coord_x;
    degree[0] = 3;
  }
  else if (dim == 2)
  {
    coord[0] = coord_x;
    coord[1] = coord_y;

    degree[0] = 3;
    degree[1] = 2;
  }
  else if (dim == 3)
  {
    coord[0] = coord_x;
    coord[1] = coord_y;
    coord[2] = coord_z;

    degree[0] = 3;
    degree[1] = 2;
    degree[2] = 1;
  }

  // Defining used types.
  static const int codim = range - dim;
  using GridType = Grid<dim>;
  using SpSpaceType = SplineSpace<dim, range, rank>;
  using RefBasisType = ReferenceBasis<dim, range, rank>;
  using BSplineType = BSpline<dim, range, rank>;
  using NURBSType = NURBS<dim, range, rank>;
  using ScalarSpSpaceType = SplineSpace<dim, 1, 1>;
  using ScalarBSplineType = BSpline<dim, 1, 1>;
  using ScalarRefBasisType = ReferenceBasis<dim, 1, 1>;
  using WeightFuncType = IgGridFunction<dim, 1>;
  using ScalarGridFuncType = GridFunction<dim, 1>;
  using GridFuncType = GridFunction<dim, range>;
  using DomainType = Domain<dim, codim>;
  using ConstGridFunc = grid_functions::ConstantGridFunction<dim, range>;
  using ConstFuncType = functions::ConstantFunction<dim, codim, range, rank>;
  using FuncType = Function<dim, codim, range, rank>;
  using PhysBasisType = PhysicalBasis<dim, range, rank, codim>;

  auto grid = GridType::create(coord);
  container->insert_const_object<GridType>(grid);

  auto ssp = SpSpaceType::create(degree,grid);
  container->insert_const_object<SpSpaceType>(ssp);

  auto  bsp = BSplineType::create(ssp);
  container->insert_const_object<RefBasisType>(bsp);

  auto scalar_basis = ScalarBSplineType::create(ScalarSpSpaceType::create(degree, grid));

  container->insert_const_object<ScalarRefBasisType>(scalar_basis);

  const auto n_scalar_basis = scalar_basis->get_num_basis();

  IgCoefficients weights;
  for (int dof = 0 ; dof < n_scalar_basis ; ++dof)
    weights[dof] = 1.0;

  const auto w_func = WeightFuncType::create(scalar_basis,weights);
  w_func->set_name("my_weight_function");

  container->insert_const_object<ScalarGridFuncType>(w_func);

  auto nurbs_basis = NURBSType::create(bsp, w_func);
  container->insert_const_object<RefBasisType>(nurbs_basis);

  Values<dim, range, 1> val;
  const auto const_grid_func = ConstGridFunc::create(grid, val);
  container->insert_const_object<GridFuncType>(const_grid_func);
  const_grid_func->set_name("my_const_grid_function");

  const auto domain = DomainType::create(const_grid_func);
  domain->set_name("my_domain");
  container->insert_const_object<DomainType>(domain);

  const auto phys_basis = PhysBasisType::create(nurbs_basis, domain);
  container->insert_const_object<PhysBasisType>(phys_basis);

  const auto const_func = ConstFuncType::create(domain, val);
  container->insert_const_object<FuncType>(const_func);
  const_func->set_name("my_const_function");
}


int main()
{

  const auto container = ObjectsContainer::create();

  insert_objects<1, 1>(container);
  insert_objects<1, 2>(container);
  insert_objects<1, 3>(container);
  insert_objects<2, 2>(container);
  insert_objects<2, 3>(container);
  insert_objects<3, 3>(container);

  OUTSTART
  container->print_info(out);
  OUTEND

  return 0;
}
