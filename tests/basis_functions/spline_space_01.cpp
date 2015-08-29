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
 *  @brief Construction of different spline spaces
 *  @author pauletti
 *  @date 2015-08-10
 */

// TODO (pauletti, Dec 26, 2014): make this test dim independent
#include "../tests.h"
#include <igatools/basis_functions/spline_space.h>


template <int dim,int range>
void print_boundary_and_repeated_knots(
  std::shared_ptr<SplineSpace<dim,range>> sp_spec,
  typename SplineSpace<dim,range>::BoundaryKnotsTable bdry_knots)
{
  using SplineSpace = SplineSpace<dim,range>;

  typename SplineSpace::EndBehaviour eb(BasisEndBehaviour::end_knots);
  typename SplineSpace::EndBehaviourTable ebt(eb);

  auto rep_knots = sp_spec->compute_knots_with_repetition(ebt, bdry_knots);
  out << "Boundary knots:\n";
  for (const auto &v : bdry_knots)
    for (const auto &w : v)
      w.print_info(out);
  out << "Repeated knots:\n";
  for (const auto &v : rep_knots)
    v.print_info(out);
}



void test_1d()
{
  OUTSTART

  const int dim=1;
  using SplineSpace = SplineSpace<dim>;
  using MultiplicityTable = typename SplineSpace::MultiplicityTable;

  auto grid = CartesianGrid<dim>::create(4);
  typename SplineSpace::DegreeTable deg {{2}};
  auto int_mult = MultiplicityTable({ {{1,3}} });
  auto sp_spec = SplineSpace::create(deg, grid, int_mult);

  CartesianProductArray<Real,2> bn_x {{-0.5, 0, 0}, {1.1, 1.2, 1.3}};
  typename SplineSpace::BoundaryKnotsTable bdry_knots { {bn_x} };

  print_boundary_and_repeated_knots(sp_spec,bdry_knots);

  OUTEND
}



void test_2d()
{
  OUTSTART

  const int dim=2;
  using SplineSpace = SplineSpace<dim>;
  using MultiplicityTable = typename SplineSpace::MultiplicityTable;
  auto grid = CartesianGrid<dim>::create({3,5});
  typename SplineSpace::DegreeTable deg {{1,3}};

  auto int_mult = MultiplicityTable({ {{1}, {1,3,1}} });

  auto sp_spec = SplineSpace::create(deg, grid, int_mult);

  iga::CartesianProductArray<double, 2> bk_x {{-0.5, 0}, {1.2, 1.3}};
  iga::CartesianProductArray<double, 2> bk_y {{-0.6,0,0,0}, {1,1.1,1.6, 1.6}};
  typename SplineSpace::BoundaryKnotsTable bdry_knots { {bk_x, bk_y} };

  print_boundary_and_repeated_knots(sp_spec,bdry_knots);

  OUTEND
}


void test_3d()
{
  OUTSTART

  const int dim=3;
  using SplineSpace = SplineSpace<dim>;
  using MultiplicityTable = typename SplineSpace::MultiplicityTable;
  auto grid = CartesianGrid<dim>::create({3,4,5});
  typename SplineSpace::DegreeTable deg {{1,3,0}};
  auto int_mult = MultiplicityTable({ {{1}, {1,3}, {1,1,1}} });

  auto sp_spec = SplineSpace::create(deg, grid, int_mult);

  iga::CartesianProductArray<double, 2> bk_x {{-0.5, 0}, {1.2, 1.3}};
  iga::CartesianProductArray<double, 2> bk_y {{-0.6,0,0,0}, {1,1,1.6, 1.6}};
  iga::CartesianProductArray<double, 2> bk_z {{-0.6}, {1.6}};
  typename SplineSpace::BoundaryKnotsTable bdry_knots { {bk_x, bk_y, bk_z} };

  print_boundary_and_repeated_knots(sp_spec,bdry_knots);

  OUTEND
}


void test_2d_2()
{
  OUTSTART

  const int dim=2;
  const int range=2;
  using SplineSpace = SplineSpace<dim, range, 1>;
  using MultiplicityTable = typename SplineSpace::MultiplicityTable;
  auto grid = CartesianGrid<dim>::create({3,4});
  typename SplineSpace::DegreeTable deg {{1,3},{3,1}};

  auto int_mult = MultiplicityTable({ {{1}, {1,3}},{{1}, {1,1}}});

  auto sp_spec = SplineSpace::create(deg, grid, int_mult);

  iga::CartesianProductArray<double, 2> bk_x {{-0.5, 0}, {1.2, 1.3}};
  iga::CartesianProductArray<double, 2> bk_y {{-0.6,0,0,0}, {1,1,1.6, 1.6}};

  typename SplineSpace::BoundaryKnotsTable bdry_knots { {bk_x, bk_y}, {bk_y, bk_x} };

  print_boundary_and_repeated_knots(sp_spec,bdry_knots);

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
