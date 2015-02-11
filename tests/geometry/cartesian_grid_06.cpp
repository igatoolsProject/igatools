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
 *  Test for CartesianGrid::find_elements_of_points
 *
 *  author: pauletti
 *  date: 2014-08-07
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>

template<int dim>
void do_test()
{
    TensorSize<dim> n_knots;
    for (int i = 0; i < dim; ++i)
        n_knots[i] = 2*i+3;
    auto grid = CartesianGrid<dim>::create(n_knots);

    out << "Dimension: " << dim << endl;
    const auto n_elems = grid->get_num_active_elems();
    for (int i = 0; i < n_elems; ++i)
    {
        const auto ti = grid->flat_to_tensor(i);
        const auto fi = grid->tensor_to_flat(ti);
        out << "Element " << fi << ": " << ti << endl;
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
