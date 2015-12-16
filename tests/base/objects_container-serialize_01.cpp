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
 *  @brief  Test for Objects container serialization
 *  @author P. Antolin
 *  @date 2015
 */

#include "../tests.h"

#include <igatools/base/objects_container.h>
#include <igatools/geometry/grid.h>

#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/nurbs_element.h>
#include <igatools/functions/ig_function.h>

void serialize_deserialize(const ObjectsContainer &container)
{
  out.begin_item("Original ObjectsContainer.");
  container.print_info(out);
  out.end_item();


  // serialize the ObjectsContainer to a xml file
  std::string filename = "objects_container.xml";
  // serialize the ObjectsContainer to an xml file
  {
    std::ofstream xml_ostream(filename);
    OArchive xml_out(xml_ostream);
    xml_out << container;
  }

  // de-serialize the ObjectsContainer from a xml file
  ObjectsContainer container_new;
  {
    std::ifstream xml_istream(filename);
    IArchive xml_in(xml_istream);
    xml_in >> container_new;
  }

  out.begin_item("ObjectsContainer after serialize-deserialize.");
  container_new.print_info(out);
  out.end_item();
  //*/
}

template< int dim>
void insert_objects(const std::shared_ptr<ObjectsContainer> container)
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
  static const int range = dim;
  static const int codim = range - dim;
  static const int rank = 1;
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
  using IgGridFunc = IgGridFunction<dim, dim+codim>;
  using IgFuncType = IgFunction<dim, codim, range, rank>;
  using FuncType = Function<dim, codim, range, rank>;
  using PhysSpaceType = PhysicalSpaceBasis<dim, range, rank, codim>;

  auto grid = GridType::create(coord);
  container->insert_const_object<GridType>(grid);

  auto ssp = SpSpaceType::create(degree,grid);
  container->insert_const_object<SpSpaceType>(ssp);

  auto  bsp = BSplineType::create(ssp);
  container->insert_const_object<RefSpaceType>(bsp);

  auto scalar_space = ScalarBSplineType::create(ScalarSpSpaceType::create(degree, grid));

  container->insert_const_object<ScalarRefSpaceType>(scalar_space);

  const auto n_scalar_basis = scalar_space->get_num_basis();

  IgCoefficients weights;
  for (int dof = 0 ; dof < n_scalar_basis ; ++dof)
    weights[dof] = 1.0;

  const auto w_func = WeightFuncType::create(scalar_space,weights);

  container->insert_const_object<ScalarGridFuncType>(w_func);

  auto nurbs_space = NURBSType::create(bsp, w_func);
  container->insert_const_object<RefSpaceType>(nurbs_space);

  Epetra_SerialComm comm;
  auto map = EpetraTools::create_map(*nurbs_space, "active", comm);
  const auto pts = EpetraTools::create_vector(*map);
  (*pts)[0] = 1.;
  auto ig_grid_func = IgGridFunc::create(nurbs_space, *pts, "active");
  container->insert_const_object<GridFuncType>(ig_grid_func);

  const auto domain = DomainType::create(ig_grid_func);
  domain->set_name("my_domain");
  container->insert_const_object<DomainType>(domain);

  const auto phys_space = PhysSpaceType::create(nurbs_space, domain);
  container->insert_const_object<PhysSpaceType>(phys_space);

  auto map_2 = EpetraTools::create_map(*phys_space, "active", comm);
  auto coeff = EpetraTools::create_vector(*map_2);
  (*coeff)[0] = 2.;
  auto ig_func = IgFuncType::create(phys_space, *coeff);
  ig_func->set_name("my_function");
  container->insert_const_object<FuncType>(ig_func);


}


int main()
{

  const auto container = ObjectsContainer::create();

  insert_objects<1>(container);
  insert_objects<2>(container);
  insert_objects<3>(container);

  OUTSTART
  serialize_deserialize(*container);
  OUTEND

  return 0;
}
