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
 *  Test for the index of a GridElement
 *
 *  author: pauletti
 *  date: 2014-08-07
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/grid_element.h>

template<int dim>
void do_test()
{
    TensorSize<dim> n_knots;
    for (int i = 0; i < dim; ++i)
        n_knots[i] = 2*i+3;
    auto grid = CartesianGrid<dim>::create(n_knots);

    out << "Dimension: " << dim << endl;

    int i = 0;
    for (const auto &elem : *grid)
    {
        out << "Element " << i << ": " << elem.get_index() << endl;
        ++i;
    }
    out << endl;
    out << endl;
}

int main()
{
    do_test<1>();
    do_test<2>();
    do_test<3>();

    return 0;
}
