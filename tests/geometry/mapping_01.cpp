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
 *  Test for the linear mapping class.
 *  This test the evaluation part.
 *  author: pauletti
 *  date: 2012-12-19
 *
 */

#include "../tests.h"

#include <igatools/geometry/mapping_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_element_accessor.h>

template <int dim, int codim >
void test_evaluate()
{
    const int space_dim = dim + codim;
    Derivatives<dim,space_dim,1,1> A;
    Points<space_dim> b;

    //Dilation
    for (int i=0; i<dim; ++i)
        A[i][i] = i+1;
    //Traslation
    for (int i=0; i<dim; ++i)
        b[i] = i+1;

    auto map = LinearMapping<dim,codim>::create(A,b);

    out << "Linear mapping" << "<" << dim << "," << space_dim << ">" << endl;
    out << "A =" << endl << A << endl;
    out << "b =" << b << endl << endl;

    QTrapez<dim> quad;

    auto elem = map->begin();
    elem->init_cache(ValueFlags::point|ValueFlags::map_gradient|ValueFlags::map_hessian, quad);
    elem->fill_cache();
    auto values=elem->get_map_values();
    auto gradients = elem->get_map_gradients();
    auto hessians = elem->get_map_hessians();

    out << "Values:" << endl;
    values.print_info(out);
    out << "Gradients:" << endl;
    gradients.print_info(out);
    out << "Hessians:" << endl;
    hessians.print_info(out);

}

int main()
{
    out.depth_console(10);

    test_evaluate<2,0>();
    test_evaluate<3,0>();
    test_evaluate<2,1>();
    test_evaluate<1,1>();
    test_evaluate<1,2>();

}
