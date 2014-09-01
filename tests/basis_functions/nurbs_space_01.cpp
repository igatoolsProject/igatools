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

#include "../tests.h"

#include <igatools/basis_functions/nurbs_space.h>


template< int dim, int range, int rank = 1>
void do_test()
{
    using iga::vector;
    vector<Real> coord_x {0,1,2,3,4};
    vector<Real> coord_y {5,6,7,8};
    vector<Real> coord_z {9, 10, 11};

    CartesianProductArray<Real, dim> coord;
    CartesianProductArray<Index , dim>  mult;
    TensorIndex<dim> degree;

    if (dim == 1)
    {
        coord.copy_data_direction(0,coord_x);
        degree[0] = 3;
    }
    else if (dim == 2)
    {
        coord.copy_data_direction(0,coord_x);
        coord.copy_data_direction(1,coord_y);

        degree[0] = 3;
        degree[1] = 2;
    }
    else if (dim == 3)
    {
        coord.copy_data_direction(0,coord_x);
        coord.copy_data_direction(1,coord_y);
        coord.copy_data_direction(2,coord_z);

        degree[0] = 3;
        degree[1] = 2;
        degree[2] = 1;
    }




    using Space = NURBSSpace< dim, range, rank >;
    using WeightsTable = typename Space::WeightsTable;
    using DegreeTable = typename Space::DegreeTable;
    auto  knots = CartesianGrid<dim>::create(coord);
    DegreeTable deg(degree);

    auto  bsp = BSplineSpace<dim, range, rank >::create(deg, knots);
    WeightsTable weights;
    const auto n_basis = bsp->get_num_basis_table();

    for (auto comp : Space::components)
        weights(comp).resize(n_basis(comp),1.0);

    for (auto &w : weights)
        w.print_info(out);

    auto nurbs_space = Space::create(deg, knots, weights);
    nurbs_space->print_info(out);
}


int main()
{
    do_test<1, 1>();
    do_test<1, 2>();
    do_test<1, 3>();

    do_test<2, 1>();
    do_test<2, 2>();
    do_test<2, 3>();

    do_test<3, 1>();
    do_test<3, 3>();

    return 0;
}
