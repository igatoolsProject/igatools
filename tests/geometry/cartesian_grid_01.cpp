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
 *  Test for CartesianGrid constructors
 *
 *  author: pauletti
 *  date: 2013-10-25
 *  QA: v0.2 (2013-10-25)
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>

#include <igatools/base/quadrature.h>

template<int dim>
void def_const()
{
    auto grid = CartesianGrid<dim>::create();
    grid->print_info(out);
    out << endl;
}

template<int dim>
void uniform_const(const int n_knots)
{
    auto grid = CartesianGrid<dim>::create(n_knots);
    grid->print_info(out);
    out << endl;
}

template<int dim>
void dim_uniform_const()
{
    TensorSize<dim> n_knots;
    for (int i = 0; i < dim; ++i)
        n_knots[i] = 2*i+2;
    auto grid = CartesianGrid<dim>::create(n_knots);
    grid->print_info(out);
    out << endl;
}


template<int dim>
void non_uniform_const()
{
    TensorSize<dim> n_knots;
    for (int i = 0; i < dim; ++i)
        n_knots[i] = 2*i+2;
    int k = 0;
    CartesianProductArray<Real, dim> knots(n_knots);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j <n_knots[i] ; ++j)
            knots.entry(i,j) = k++;

    auto grid = CartesianGrid<dim>::create(knots);
    grid->print_info(out);
    out << endl;
}

int main()
{
    def_const<0>();
    def_const<1>();
    def_const<2>();
    def_const<3>();

    uniform_const<0>(3);
    uniform_const<1>(3);
    uniform_const<2>(3);
    uniform_const<3>(3);

    dim_uniform_const<0>();
    dim_uniform_const<1>();
    dim_uniform_const<2>();
    dim_uniform_const<3>();

    non_uniform_const<0>();
    non_uniform_const<1>();
    non_uniform_const<2>();
    non_uniform_const<3>();

    return 0;
}
