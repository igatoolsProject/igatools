//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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
 *  @brief  serialization of the Splinespace
 *  @author martinelli
 *  @date 2015-08-10
 */

// TODO (pauletti, Dec 26, 2014): make this test dim independent

#include "spline_space_tests_common.h"



template <int dim,int range>
void serialize_deserialize(std::shared_ptr<SplineSpace<dim,range>> spline_space)
{
  using std::to_string;
  const std::string class_name = "SplineSpace<" + to_string(dim) + "," + to_string(range) + ",1>";
  out.begin_item("Original " + class_name);
  spline_space->print_info(out);
  out.end_item();



  std::string filename = "spline_space_dim" + std::to_string(dim) + "_range" + std::to_string(range) + ".xml";
  std::string tag_name = "SplineSpace_dim" + std::to_string(dim) + "_range" + std::to_string(range);
  {
    // serialize the SplineSpace object to an xml file
    std::ofstream xml_ostream(filename);
    OArchive xml_out(xml_ostream);

    xml_out << spline_space;
  }

  auto grid = Grid<dim>::create();
  typename SplineSpace<dim,range>::DegreeTable deg(TensorIndex<dim>(3));
  typename SplineSpace<dim,range>::MultiplicityTable mult;
  auto spline_space_new = SplineSpace<dim,range>::create(deg, grid, mult);
  {
    // de-serialize the SplineSpace object from an xml file
    std::ifstream xml_istream(filename);
    IArchive xml_in(xml_istream);
    xml_in >> spline_space_new;
  }
  out.begin_item(class_name + " after serialize-deserialize.");
  spline_space_new->print_info(out);
  out.end_item();
  //*/
}


void test_1d()
{
  OUTSTART

  const int dim=1;
  using SplineSpace = SplineSpace<dim>;
  using MultiplicityTable = typename SplineSpace::MultiplicityTable;

  auto grid = Grid<dim>::create(4);
  typename SplineSpace::DegreeTable deg {{2}};
  auto int_mult = MultiplicityTable({ {{1,3}} });
  auto spline_space = SplineSpace::create(deg, grid, int_mult);

  serialize_deserialize(spline_space);

  SafeSTLArray<SafeSTLVector<Real>,2> bn_x {{-0.5, 0, 0}, {1.1, 1.2, 1.3}};
  typename SplineSpace::BoundaryKnotsTable bdry_knots { {bn_x} };

  typename SplineSpace::EndBehaviour eb(BasisEndBehaviour::end_knots);
  typename SplineSpace::EndBehaviourTable ebt(eb);

  print_boundary_and_repeated_knots(*spline_space,bdry_knots,ebt);

  OUTEND
}



void test_2d()
{
  OUTSTART

  const int dim=2;
  using SplineSpace = SplineSpace<dim>;
  using MultiplicityTable = typename SplineSpace::MultiplicityTable;
  auto grid = Grid<dim>::create({3,5});
  typename SplineSpace::DegreeTable deg {{1,3}};

  auto int_mult = MultiplicityTable({ {{1}, {1,3,1}} });

  auto spline_space = SplineSpace::create(deg, grid, int_mult);

  serialize_deserialize(spline_space);

  SafeSTLArray<SafeSTLVector<Real>,2> bk_x {{-0.5, 0}, {1.2, 1.3}};
  SafeSTLArray<SafeSTLVector<Real>,2> bk_y {{-0.6,0,0,0}, {1,1.1,1.6, 1.6}};
  typename SplineSpace::BoundaryKnotsTable bdry_knots { {bk_x, bk_y} };

  typename SplineSpace::EndBehaviour eb(BasisEndBehaviour::end_knots);
  typename SplineSpace::EndBehaviourTable ebt(eb);

  print_boundary_and_repeated_knots(*spline_space,bdry_knots,ebt);

  OUTEND
}


void test_3d()
{
  OUTSTART

  const int dim=3;
  using SplineSpace = SplineSpace<dim>;
  using MultiplicityTable = typename SplineSpace::MultiplicityTable;
  auto grid = Grid<dim>::create({3,4,5});
  typename SplineSpace::DegreeTable deg {{1,3,0}};
  auto int_mult = MultiplicityTable({ {{1}, {1,3}, {1,1,1}} });

  auto spline_space = SplineSpace::create(deg, grid, int_mult);

  serialize_deserialize(spline_space);

  SafeSTLArray<SafeSTLVector<Real>,2> bk_x {{-0.5, 0}, {1.2, 1.3}};
  SafeSTLArray<SafeSTLVector<Real>,2> bk_y {{-0.6,0,0,0}, {1,1,1.6, 1.6}};
  SafeSTLArray<SafeSTLVector<Real>,2> bk_z {{-0.6}, {1.6}};
  typename SplineSpace::BoundaryKnotsTable bdry_knots { {bk_x, bk_y, bk_z} };

  typename SplineSpace::EndBehaviour eb(BasisEndBehaviour::end_knots);
  typename SplineSpace::EndBehaviourTable ebt(eb);

  print_boundary_and_repeated_knots(*spline_space,bdry_knots,ebt);

  OUTEND
}


void test_2d_2()
{
  OUTSTART

  const int dim=2;
  const int range=2;
  using SplineSpace = SplineSpace<dim, range, 1>;
  using MultiplicityTable = typename SplineSpace::MultiplicityTable;
  auto grid = Grid<dim>::create({3,4});
  typename SplineSpace::DegreeTable deg {{1,3},{3,1}};

  auto int_mult = MultiplicityTable({ {{1}, {1,3}},{{1}, {1,1}}});

  auto spline_space = SplineSpace::create(deg, grid, int_mult);

  serialize_deserialize(spline_space);

  SafeSTLArray<SafeSTLVector<Real>,2> bk_x {{-0.5, 0}, {1.2, 1.3}};
  SafeSTLArray<SafeSTLVector<Real>,2> bk_y {{-0.6,0,0,0}, {1,1,1.6, 1.6}};

  typename SplineSpace::BoundaryKnotsTable bdry_knots { {bk_x, bk_y}, {bk_y, bk_x} };

  typename SplineSpace::EndBehaviour eb(BasisEndBehaviour::end_knots);
  typename SplineSpace::EndBehaviourTable ebt(eb);

  print_boundary_and_repeated_knots(*spline_space,bdry_knots,ebt);

  OUTEND
}


int main()
{
  out.depth_console(10);

  test_1d();
  test_2d();
  test_3d();

  test_2d_2();

  return 0;
}
