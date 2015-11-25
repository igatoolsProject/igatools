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
 *  Test for  Gauss-Legendre quadrature scheme
 *
 *  author: pauletti
 *  date: 2015-03-06
 *
 */

#include "../tests.h"

#include <igatools/base/objects_container.h>
#include <igatools/geometry/grid.h>

#include <igatools/geometry/grid_function_lib.h>
#include <igatools/functions/function_lib.h>
#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/nurbs_element.h>

template< int dim, int range, int rank = 1>
void insert_objects (const std::shared_ptr<ObjectsContainer> container,
                     Index &object_id)
{
  using iga::SafeSTLVector;
  SafeSTLVector<Real> coord_x {0,1,2,3,4};
  SafeSTLVector<Real> coord_y {5,6,7,8};
  SafeSTLVector<Real> coord_z {9, 10, 11};

  SafeSTLArray<SafeSTLVector<Real>, dim> coord;
  CartesianProductArray<Index , dim>  mult;
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
  using RefSpaceType = ReferenceSpaceBasis<dim, range, rank>;
  using BSplineType = BSpline<dim, range, rank>;
  using NURBSType = NURBS<dim, range, rank>;
  using ScalarSpSpaceType = SplineSpace<dim, 1, 1>;
  using ScalarBSplineType = BSpline<dim, 1, 1>;
  using ScalarRefSpaceType = ReferenceSpaceBasis<dim, 1, 1>;
  using WeightFuncType = IgGridFunction<dim, 1>;
  using ScalarGridFuncType = GridFunction<dim, 1>;
  using GridFuncType = GridFunction<dim, range>;
  using DomainType = Domain<dim, codim>;
  using ConstGridFunc = grid_functions::ConstantGridFunction<dim, range>;
  using ConstFuncType = functions::ConstantFunction<dim, codim, range, rank>;
  using FuncType = Function<dim, codim, range, rank>;
  using PhysSpaceType = PhysicalSpaceBasis<dim, range, rank, codim>;

  auto grid = GridType::create(coord);
  container->insert_object<GridType>(grid, object_id);
  ++object_id;

  auto ssp = SpSpaceType::create(degree,grid);
  container->insert_object<SpSpaceType>(ssp, object_id);
  ++object_id;

  auto  bsp = BSplineType::create(ssp);
  container->insert_object<RefSpaceType>(bsp, object_id);
  ++object_id;

  auto scalar_space = ScalarBSplineType::create(ScalarSpSpaceType::create(degree,grid));

  container->insert_object<ScalarRefSpaceType>(scalar_space, object_id);
  ++object_id;

  const auto n_scalar_basis = scalar_space->get_num_basis();

  IgCoefficients weights;
  for (int dof = 0 ; dof < n_scalar_basis ; ++dof)
    weights[dof] = 1.0;

  const auto w_func = WeightFuncType::create(scalar_space,weights);

  container->insert_object<ScalarGridFuncType>(w_func, object_id);
  ++object_id;

  auto nurbs_space = NURBSType::create(bsp, w_func);
  container->insert_object<RefSpaceType>(nurbs_space, object_id);
  ++object_id;

  const Values<dim, range, 1> val;
  const auto const_grid_func = ConstGridFunc::create(grid, val);
  container->insert_object<GridFuncType>(const_grid_func, object_id);
  ++object_id;

  const auto domain = DomainType::create (const_grid_func, "my_domain");
  container->insert_object<DomainType>(domain, object_id);
  ++object_id;

  const auto phys_space = PhysSpaceType::create(nurbs_space, domain);
  container->insert_object<PhysSpaceType>(phys_space, object_id);
  ++object_id;

  const auto const_func = ConstFuncType::create(domain, val);
  container->insert_object<FuncType>(const_func, object_id);
  ++object_id;
}


template< int dim, int range, int rank = 1>
void retrieve_objects(const std::shared_ptr<ObjectsContainer> container,
                      Index &object_id)
{
  OUTSTART

  // Defining used types.
  static const int codim = range - dim;
  using GridType = Grid<dim>;
  using SpSpaceType = SplineSpace<dim, range, rank>;
  using RefSpaceType = ReferenceSpaceBasis<dim, range, rank>;
  using BSplineType = BSpline<dim, range, rank>;
  using NURBSType = NURBS<dim, range, rank>;
  using ScalarBSplineType = BSpline<dim, 1, 1>;
  using ScalarRefSpaceType = ReferenceSpaceBasis<dim, 1, 1>;
  using ScalarGridFuncType = GridFunction<dim, 1>;
  using GridFuncType = GridFunction<dim, range>;
  using DomainType = Domain<dim, codim>;
  using FuncType = Function<dim, codim, range, rank>;
  using PhysSpaceType = PhysicalSpaceBasis<dim, range, rank, codim>;

  const auto grid = container->get_object<GridType>(object_id);
  ++object_id;
  grid->print_info(out);
  out << endl;

  const auto ssp = container->get_object<SpSpaceType>(object_id);
  ++object_id;
  ssp->print_info(out);
  out << endl;

  const auto rs0 = container->get_object<RefSpaceType>(object_id);
  const auto bsp0 = std::dynamic_pointer_cast<BSplineType>(rs0);
  ++object_id;
  rs0->print_info(out);
  out << endl;
  bsp0->print_info(out);
  out << endl;

  const auto rs1 = container->get_object<ScalarRefSpaceType>(object_id);
  const auto bsp1 = std::dynamic_pointer_cast<ScalarBSplineType>(rs1);
  ++object_id;
  rs1->print_info(out);
  out << endl;
  bsp1->print_info(out);
  out << endl;

  const auto gf0 = container->get_object<ScalarGridFuncType>(object_id);
  ++object_id;
  gf0->print_info(out);
  out << endl;

  const auto rs2 = container->get_object<RefSpaceType>(object_id);
  const auto nr = std::dynamic_pointer_cast<NURBSType>(rs2);
  ++object_id;
  rs2->print_info(out);
  out << endl;
  nr->print_info(out);
  out << endl;

  const auto gf1 = container->get_object<GridFuncType>(object_id);
  ++object_id;
  gf1->print_info(out);
  out << endl;

  const auto dm = container->get_object<DomainType>(object_id);
  ++object_id;
  dm->print_info(out);
  out << endl;

  const auto ps = container->get_object<PhysSpaceType>(object_id);
  ++object_id;
  ps->print_info(out);
  out << endl;

  const auto fn = container->get_object<FuncType>(object_id);
  ++object_id;
  //  fn->print_info(out);
  out << endl;

  OUTEND

}


int main()
{

  const auto container = ObjectsContainer::create();

  Index object_id = 0;
  insert_objects<1, 1>(container, object_id);
  insert_objects<1, 2>(container, object_id);
  insert_objects<1, 3>(container, object_id);
  insert_objects<2, 2>(container, object_id);
  insert_objects<2, 3>(container, object_id);
  insert_objects<3, 3>(container, object_id);

  object_id = 0;
  retrieve_objects<1, 1>(container, object_id);
  retrieve_objects<1, 2>(container, object_id);
  retrieve_objects<1, 3>(container, object_id);
  retrieve_objects<2, 2>(container, object_id);
  retrieve_objects<2, 3>(container, object_id);
  retrieve_objects<3, 3>(container, object_id);

  return 0;
}
