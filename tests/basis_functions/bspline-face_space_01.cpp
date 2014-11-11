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
 *  Test for the BsplineSpace class face subspace extraction function.
 *  Here we print the information of the face spaces thus extracted.
 *  author: pauletti
 *  date: Jan 29, 2013
 *  updated: 2013-04-02
 */

#include "../tests.h"
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/new_bspline_space.h>
//#include <igatools/basis_functions/space_tools.h>

template<int k, int dim, int range=1, int rank=1>
void run_test(TensorSize<dim> n, const int degree = 1)
{
    using Space = NewBSplineSpace<dim, range, rank>;

    auto grid = CartesianGrid<dim>::create(n);
    auto space = Space::create(degree, grid);

    typename Space::template InterSpaceMap<k> dof_map;
    typename CartesianGrid<dim>::template InterGridMap<k> elem_map;

    for (auto i : UnitElement<dim>::template elems_ids<k>())
    {

        out << "Face space: " << i << endl;
        auto sub_space =
            space->template get_ref_sub_space<k>(i, dof_map);
        sub_space->print_info(out);

        out << "Dofs face to space mapping: " << endl;
        dof_map.print_info(out);
        out << endl;
    }
}



int main()
{
    //out.depth_console(10);

    run_test<0,1>(TensorSize<1>(arr::sequence<1>(2)));
    run_test<1,2>(TensorSize<2>(arr::sequence<2>(2)));
    //run_test<1,2,3>(TensorSize<2>(arr::sequence<2>(2)));
    run_test<2,3>(TensorSize<3>(arr::sequence<3>(2)));

    return  0;
}
