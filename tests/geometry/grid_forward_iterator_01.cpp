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
 *  Test iterators and accessors on a CartesianGrid
 *  pauletti 2012-11-18
 *
 */


#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>

#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/geometry/grid_forward_iterator.h>
#include <igatools/geometry/unit_element.h>


using namespace std ;


// Default constructor
template <int dim>
void test_accessor1()
{
    CartesianGrid< dim > cartesian_grid ;
    GridForwardIterator<CartesianGridElementAccessor<dim> > element(cartesian_grid, 0);

    for (int i = 0; i < UnitElement<dim>::vertices_per_element; ++i)
    {
        out << element->vertex(i) << endl;
    }
}


// Default constructor using begin and end
template <int dim>
void test_accessor2()
{
    CartesianGrid< dim > cartesian_grid ;



    typename CartesianGrid<dim>::ElementIterator element = cartesian_grid.begin();
    typename CartesianGrid<dim>::ElementIterator end_element = cartesian_grid.end();
    for (; element != end_element; ++element)
        for (int i = 0; i < UnitElement<dim>::vertices_per_element; ++i)
            out<<element->vertex(i) << endl;
}




template <int dim>
void test_accessor3()
{
    CartesianProductArray< Real, dim > coords;

    const int n_knots = 4;

    for (int i = 0; i < dim; ++i)
    {
        vector<Real> tmp_coord;
        for (int j = 0; j < n_knots; ++j)
            tmp_coord.push_back(j);
        coords.copy_data_direction(i,tmp_coord);
    }
    CartesianGrid< dim > cartesian_grid(coords) ;

    GridForwardIterator<CartesianGridElementAccessor<dim> > element(cartesian_grid, 0);

    const int n_elements = cartesian_grid.get_num_elements();
    for (int j = 0; j < n_elements; ++j, ++element)
    {
        out << "Element: " << j << endl;
        for (int i = 0; i < UnitElement<dim>::vertices_per_element; ++i)
        {
            out << "\t" << element->vertex(i) << endl;

        }
    }


}




template <int dim>
void test_accessor4()
{
    CartesianProductArray< Real, dim> coords;

    const int n_knots = 4;

    for (int i = 0; i < dim; ++i)
    {
        vector<Real> tmp_coord;
        for (int j = 0; j < n_knots; ++j)
            tmp_coord.push_back(j);
        coords.copy_data_direction(i,tmp_coord);
    }


    CartesianGrid< dim > cartesian_grid(coords) ;

    typename CartesianGrid<dim>::ElementIterator element = cartesian_grid.begin();
    typename CartesianGrid<dim>::ElementIterator end_element = cartesian_grid.end();
    for (; element != end_element; ++element)
    {
        //element->indices();
        for (int i = 0; i < UnitElement<dim>::vertices_per_element; ++i)
            out<<element->vertex(i) << endl;
    }
}



int main(int argc, char *argv[])
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
