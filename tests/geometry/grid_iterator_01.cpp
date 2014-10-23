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
 *  Test the CartesianGrid iterator accessor .
 *  author: pauletti
 *  date:  2014-04-24
 *
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/geometry/grid_forward_iterator.h>
#include <igatools/geometry/unit_element.h>


// Accessor constructor
template <int dim>
void test_accessor1()
{
    auto grid = CartesianGrid<dim>::create();
    GridForwardIterator<CartesianGridElement<dim> > elem(grid, 0);

    for (int i = 0; i < UnitElement<dim>::vertices_per_element; ++i)
        out << elem->vertex(i) << endl;
}



// Iterator for default grid
template <int dim>
void test_accessor2()
{
    auto grid    = CartesianGrid<dim>::create();
    auto elem = grid->begin();
    auto end     = grid->end();
    for (; elem != end; ++elem)
        for (int i = 0; i < UnitElement<dim>::vertices_per_element; ++i)
            out<<elem->vertex(i) << endl;
}



// Iterator for user coords grid
template <int dim>
void test_accessor3()
{
    CartesianProductArray<Real, dim> coords;
    const int n_knots = 4;
    for (int i = 0; i < dim; ++i)
    {
        iga::vector<Real> tmp_coord;
        for (int j = 0; j < n_knots; ++j)
            tmp_coord.push_back(j);
        coords.copy_data_direction(i,tmp_coord);
    }

    auto grid = CartesianGrid<dim>::create(coords);

    GridForwardIterator<CartesianGridElement<dim> > elem(grid, 0);

    const int n_elements = grid->get_num_active_elems();
    for (int j = 0; j < n_elements; ++j, ++elem)
    {
        out << "Element: " << j << endl;
        for (int i = 0; i < UnitElement<dim>::vertices_per_element; ++i)
        {
            out << "\t" << elem->vertex(i) << endl;

        }
    }
}



template <int dim>
void test_accessor4()
{

    CartesianProductArray<Real, dim> coords;
    const int n_knots = 4;
    for (int i = 0; i < dim; ++i)
    {
        iga::vector<Real> tmp_coord;
        for (int j = 0; j < n_knots; ++j)
            tmp_coord.push_back(j);
        coords.copy_data_direction(i,tmp_coord);
    }

    auto grid = CartesianGrid<dim>::create(coords);

    auto elem = grid->begin();
    auto end     = grid->end();
    for (; elem != end; ++elem)
        for (int i = 0; i < UnitElement<dim>::vertices_per_element; ++i)
            out<<elem->vertex(i) << endl;
}



int main()
{

    test_accessor1<1>();
    test_accessor1<2>();
    test_accessor1<3>();


    test_accessor2<1>();
    test_accessor2<2>();
    test_accessor2<3>();

    test_accessor3<1>();
    test_accessor3<2>();
    test_accessor3<3>();

    test_accessor4<1>();
    test_accessor4<2>();
    test_accessor4<3>();

    return 0;
}
