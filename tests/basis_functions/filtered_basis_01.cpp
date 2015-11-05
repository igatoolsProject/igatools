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
 *  Test for building a matrix on a space of an igfunction
 *
 *  author: pauletti
 *  date: 2015-03-17
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>

struct DofProp
{
  static const std::string interior;
  static const std::string dirichlet;
  static const std::string neumman;
};

const std::string DofProp::interior = "interior";
const std::string DofProp::dirichlet  = "dirichlet";
const std::string DofProp::neumman  = "neumman";

template<int dim, int range = 1, int rank = 1>
void filtered_dofs(const int deg = 1, const int n_knots = 3)
{
  OUTSTART

  using Basis = BSpline<dim, range, rank>;

  auto grid = Grid<dim>::create(n_knots);
  auto space = SplineSpace<dim,range,rank>::create(deg,grid);
  auto basis = Basis::create(space);
  auto dof_dist = space->get_dof_distribution();
  dof_dist->add_dofs_property(DofProp::interior);
  dof_dist->add_dofs_property(DofProp::dirichlet);
  dof_dist->add_dofs_property(DofProp::neumman);

  std::set<Index> int_dofs= {4};
  dof_dist->set_dof_property_status(DofProp::interior, int_dofs,true);
  std::set<Index> dir_dofs= {6,3,0, 1, 2, 5, 8};
  dof_dist->set_dof_property_status(DofProp::dirichlet, dir_dofs,true);
  std::set<Index> neu_dofs= {7};
  dof_dist->set_dof_property_status(DofProp::neumman, neu_dofs,true);

  auto elem = basis->begin();
  auto end  = basis->end();

  for (; elem != end; ++elem)
  {
    out << "Interior dofs:" << endl;
    out << "Number: " << elem->get_num_basis(DofProp::interior) << endl;
    elem->get_local_to_global(DofProp::interior).print_info(out);
    out << endl;

    out << "dirichlet dofs:" << endl;
    elem->get_local_to_global(DofProp::dirichlet).print_info(out);
    out << endl;

    out << "neumman dofs:" << endl;
    elem->get_local_to_global(DofProp::neumman).print_info(out);
    out << endl;
    out << endl;
  }

  OUTEND
}




int main()
{
  const int dim = 2;
  filtered_dofs<dim>();

  return 0;
}
