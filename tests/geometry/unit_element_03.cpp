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
 *  Test of unit element
 *  author: pauletti
 *  date: 2014-08-18
 *
 */

#include "../tests.h"

#include <igatools/geometry/unit_element.h>


template<int dim, int k>
void all_cube_elements()
{
    OUTSTART

    const auto size = UnitElement<dim>::template num_elem<k>();
    out << "Number of elements: " << size << endl;
    for (auto i=0; i<size; ++i)
    {
        out.begin_item("Element: " + std::to_string(i));
        auto &k_elem = UnitElement<dim>::template get_elem<k>(i);
        const int n_dir = k_elem.constant_directions.size();
        out << "Constant directions" << endl;
        for (int j=0; j<n_dir; ++j)
        {
            out << "x["<< k_elem.constant_directions[j]<< "]";
            out << " = " << k_elem.constant_values[j] << endl;
        }
        out << "Active directions" << endl;
        for (auto &dir : k_elem.active_directions)
            out << "x[" << dir << "]" << endl;
        out.end_item();
    }
    OUTEND
}

int main()
{
    out.depth_console(20);
    all_cube_elements<0,0>();
    all_cube_elements<1,1>();
    all_cube_elements<1,0>();
    all_cube_elements<2,2>();
    all_cube_elements<2,1>();
    all_cube_elements<2,0>();
    all_cube_elements<3,3>();
    all_cube_elements<3,2>();
    all_cube_elements<3,1>();
    all_cube_elements<3,0>();

    return 0;
}
