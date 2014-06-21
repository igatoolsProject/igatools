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
 *  Test for flat element id retrieval from point coordinate
 *
 *  author: martinelli
 *  date: 2014-04-24
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>

//#include <igatools/base/quadrature.h>



template<int dim>
void do_test()
{
    TensorSize<dim> n_knots;
    for (int i = 0; i < dim; ++i)
        n_knots(i) = 2*i+2;
    auto grid = CartesianGrid<dim>::create(n_knots);


    Points<dim> p_origin; // origin

    Points<dim> p_mid; // mid-point
    for (int i=0 ; i < dim ; ++i)
        p_mid(i) = 0.5;

    Points<dim> p_end; // end-point
    for (int i=0 ; i < dim ; ++i)
        p_end(i) = 1.0;


    grid->print_info(out);

    out << "The element containing p_origin="<< p_origin << " has the flat-id="
        << grid->get_element_flat_id_from_point(p_origin) <<endl;
    out << "The element containing p_mid="<< p_mid << " has the flat-id="
        << grid->get_element_flat_id_from_point(p_mid) <<endl;
    out << "The element containing p_end="<< p_end << " has the flat-id="
        << grid->get_element_flat_id_from_point(p_end) <<endl;
    out << endl;
}

int main()
{

    do_test<0>();
    do_test<1>();
    do_test<2>();
    do_test<3>();

    return 0;
}
