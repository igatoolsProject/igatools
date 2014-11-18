//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
 *  Test for CartesianGrid normal spaces
 *
 *  author: pauletti
 *  date: 2014-10-27
 */

#include "../tests.h"
#include <igatools/geometry/cartesian_grid.h>

template <int dim, int k>
void normal_space()
{
    OUTSTART

    auto grid = CartesianGrid<dim>::create();
    const int n_elems = UnitElement<dim>::template num_elem<k>();
    for (int j = 0; j < n_elems; ++j)
    {
        out << "Sub element index: " << j << endl;
        auto normals = grid->template get_normal_space<k>(j);
        out << "Outer normals: ";
        for (int i=0; i<dim-k; ++i)
            out << normals[i] << " " << endl;
    }

    OUTEND
}



int main()
{
    normal_space<1, 0>();
    normal_space<2, 1>();
    normal_space<3, 2>();

    return 0;
}
