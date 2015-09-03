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
 *  @brief  grid serialization
 *  @author pauletti
 *  @date   2015-08-19
 */

#include "../tests.h"
#include <igatools/geometry/grid.h>

template <int dim>
void serialize_deserialize(std::shared_ptr<const Grid<dim>> grid,
                           const std::string &filename)
{
  out.begin_item("Original grid.");
  grid->print_info(out);
  out.end_item();

  string tag_name = "Grid_" + std::to_string(dim) + "d";
  {
    std::ofstream xml_ostream(filename);
    OArchive xml_out(xml_ostream);

    xml_out << boost::serialization::make_nvp(tag_name.c_str(),*grid);
    xml_ostream.close();
  }

  auto grid_new = Grid<dim>::create(4);
  {
    ifstream xml_istream(filename);
    IArchive xml_in(xml_istream);
    xml_in >> BOOST_SERIALIZATION_NVP(*grid_new);
    xml_istream.close();
  }

  out.begin_item("Grid after serialize-deserialize.");
  grid_new->print_info(out);
  out.end_item();
}


template<int dim>
void serialize_grid(const int n_knots = 4)
{
  OUTSTART

  auto grid = Grid<dim>::const_create(n_knots);
  string filename = "grid_" + std::to_string(dim) + "d.xml";
  serialize_deserialize(grid,filename);

  OUTEND
}



int main()
{
  serialize_grid<0>();
  serialize_grid<1>();
  serialize_grid<2>();
  serialize_grid<3>();

  return 0;
}
