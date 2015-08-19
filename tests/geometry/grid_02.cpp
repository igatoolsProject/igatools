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
 *  @brief boundary ids (get and set)
 *  @author pauletti
 *  @date   2015-08-19
 */

#include "../tests.h"
#include <igatools/geometry/cartesian_grid.h>

template <int dim>
void boundary_ids()
{
    OUTSTART

    auto grid = CartesianGrid<dim>::create();

    for (auto &j : UnitElement<dim>::faces)
    {
        out << "Face number: " << j << endl;
        out << "Face boundary id: " << grid->get_boundary_id(j) << endl;
    }

    for (auto &j : UnitElement<dim>::faces)
        grid->set_boundary_id(j,j);

    for (auto &j : UnitElement<dim>::faces)
    {
        out << "Face number: " << j << endl;
        out << "Face boundary id: " << grid->get_boundary_id(j) << endl;
    }
    OUTEND
}



int main()
{
    boundary_ids<0>();
    boundary_ids<1>();
    boundary_ids<2>();
    boundary_ids<3>();

    return 0;
}
