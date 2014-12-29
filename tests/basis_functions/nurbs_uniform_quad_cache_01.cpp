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
 *  Test for the BSplineSpace UniformQuadCache
 *
 *  author: pauletti
 *  date: Aug 21, 2014
 *
 */

#include "../tests.h"

#include <igatools/basis_functions/nurbs_uniform_quad_cache.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/nurbs_element_accessor.h>




template <int dim, int range=1, int rank=1>
void uniform_space_cache(const ValueFlags flag)
{
    OUTSTART

    using iga::vector;
    vector<Real> coord_x {0,1,2,3,4};
    vector<Real> coord_y {5,6,7,8};
    vector<Real> coord_z {9, 10, 11};

    CartesianProductArray<Real, dim> coord;
    CartesianProductArray<Index,dim>  mult;
    TensorIndex<dim> deg;

    if (dim == 1)
    {
        coord.copy_data_direction(0,coord_x);
        deg[0] = 3;
    }
    else if (dim == 2)
    {
        coord.copy_data_direction(0,coord_x);
        coord.copy_data_direction(1,coord_y);

        deg[0] = 3;
        deg[1] = 2;
    }
    else if (dim == 3)
    {
        coord.copy_data_direction(0,coord_x);
        coord.copy_data_direction(1,coord_y);
        coord.copy_data_direction(2,coord_z);

        deg[0] = 3;
        deg[1] = 2;
        deg[2] = 1;
    }


    using Space = NURBSSpace< dim, range, rank >;
    using WeightsTable = typename Space::WeightsTable;
    using DegreeTable = typename Space::DegreeTable;
    auto  knots = CartesianGrid<dim>::create(coord);
    DegreeTable degree(deg);

    auto  bsp = BSplineSpace<dim, range, rank >::create(degree, knots);
    WeightsTable weights;
    const auto n_basis = bsp->get_num_basis_table();

    for (auto comp : Space::components)
        weights[comp].resize(n_basis[comp],1.0);

    auto nurbs_space = Space::create(degree, knots, weights);



    auto quad = QGauss<dim>(2);
    NURBSUniformQuadCache<dim, range, rank> cache(nurbs_space, flag, quad);
    cache.print_info(out);

    OUTEND
}



int main()
{
    out.depth_console(10);

    uniform_space_cache<1>(ValueFlags::value);
    uniform_space_cache<2>(ValueFlags::value);

    uniform_space_cache<1>(ValueFlags::gradient);
    return  0;
}
