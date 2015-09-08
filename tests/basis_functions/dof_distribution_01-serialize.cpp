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
 *  @brief  DofDistribution class serialization
 *  @author martinelli
 *  @date 2015-05-06
 *  @todo make test dimension independent and split
 */

#include "../tests.h"
#include <igatools/basis_functions/dof_distribution.h>

template <int dim>
void serialize_deserialize(const DofDistribution<dim> &dof_admin)
{
  out.begin_item("Original DofDistribution.");
  dof_admin.print_info(out);
  out.end_item();



  std::string filename = "dof_distribution_dim" + std::to_string(dim) + ".xml";
  std::string tag_name = "DofDistribution_dim" + std::to_string(dim);
  {
    // serialize the DofDistribution object to an xml file
    std::ofstream xml_ostream(filename);
    OArchive xml_out(xml_ostream);

    xml_out << dof_admin;
  }

  DofDistribution<dim> dof_admin_new;
  {
    // de-serialize the DofDistribution object from an xml file
    std::ifstream xml_istream(filename);
    IArchive xml_in(xml_istream);
    xml_in >> dof_admin_new;
  }
  out.begin_item("DofDistribution after serialize-deserialize.");
  dof_admin_new.print_info(out);
  out.end_item();
  //*/
}


template <int dim>
void dof_distribution_serialization()
{
  OUTSTART
  using SplineSpace = SplineSpace<dim>;

  typename SplineSpace::DegreeTable deg(TensorIndex<dim>(2));

  auto grid = Grid<dim>::create(4);
  auto int_mult = SplineSpace::get_multiplicity_from_regularity(InteriorReg::maximum,
                  deg, grid->get_num_intervals());
  auto sp_spec = SplineSpace::create(deg, grid,int_mult);


  auto n_basis = sp_spec->get_num_basis_table();
  auto degree = sp_spec->get_degree_table();

  DofDistribution<dim> basis_index(n_basis, degree, sp_spec->get_periodic_table());

  //-----------------------------------------------------------------
  const auto &dofs_view = basis_index.get_dofs_view();
  const std::string property_active = "active";
  // basis_index.add_dofs_property(property_active);

  for (const auto &dof : dofs_view)
  {
    if (dof % 2 == 0)
      basis_index.set_dof_property_status(property_active, dof, true);
    else
      basis_index.set_dof_property_status(property_active, dof, false);
  }
  //-----------------------------------------------------------------

  serialize_deserialize(basis_index);

  OUTEND
}

int main()
{
  out.depth_console(10);

//  dof_distribution_serialization<0>();
  dof_distribution_serialization<1>();
  dof_distribution_serialization<2>();
  dof_distribution_serialization<3>();

  return 0;
}
