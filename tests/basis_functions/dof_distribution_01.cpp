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
 *  @brief  DofDistribution class
 *  @author pauletti
 *  @date 2014-08-26
 *  @todo make test dimension independent and split
 */

#include "../tests.h"
#include <igatools/basis_functions/dof_distribution.h>


template <int dim>
void test1()
{
  OUTSTART
  using SplineSpace = SplineSpace<dim>;
  using MultiplicityTable = typename SplineSpace::MultiplicityTable;

  typename SplineSpace::DegreeTable deg {{2}};

  auto grid = Grid<dim>::const_create(4);

  auto int_mult = MultiplicityTable({ {{1,3}} });
  auto sp_spec = SplineSpace::const_create(deg, grid, int_mult);


  auto n_basis = sp_spec->get_num_basis_table();
  auto degree = sp_spec->get_degree_table();

  DofDistribution<dim> dof_admin(n_basis, degree, sp_spec->get_periodic_table());

  //-----------------------------------------------------------------
  const auto &dofs_view = dof_admin.get_dofs_view();
  const std::string property_active = "active";

  for (const auto &dof : dofs_view)
  {
    if (dof % 2 == 0)
      dof_admin.set_dof_property_status(property_active, dof, true);
    else
      dof_admin.set_dof_property_status(property_active, dof, false);
  }
  //-----------------------------------------------------------------

  dof_admin.print_info(out);

  OUTEND
}

template <int dim>
void test2()
{
  OUTSTART
  using SplineSpace = SplineSpace<dim>;

  typename SplineSpace::DegreeTable deg {{1,2}};

  auto grid = Grid<dim>::const_create({4,3});
  auto int_mult = SplineSpace::get_multiplicity_from_regularity(InteriorReg::maximum,
                  deg, grid->get_num_intervals());
  auto sp_spec = SplineSpace::const_create(deg, grid, int_mult);


  auto n_basis = sp_spec->get_num_basis_table();
  auto degree = sp_spec->get_degree_table();

  DofDistribution<dim> basis_index(n_basis, degree, sp_spec->get_periodic_table());

  //-----------------------------------------------------------------
  const auto &dofs_view = basis_index.get_dofs_view();
  const std::string property_active = "active";

  for (const auto &dof : dofs_view)
  {
    if (dof % 2 == 0)
      basis_index.set_dof_property_status(property_active, dof, true);
    else
      basis_index.set_dof_property_status(property_active, dof, false);
  }
  //-----------------------------------------------------------------

  basis_index.print_info(out);

  OUTEND
}

int main()
{
  out.depth_console(10);
  test1<1>();
  test2<2>();

  return 0;
}
