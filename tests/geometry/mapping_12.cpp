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
 *  Test for linear mapping class
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/mapping_element.h>
#include <igatools/../../source/geometry/grid_forward_iterator.cpp>
#include <igatools/base/function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_element.h>

template<int dim, int codim, int sub_dim = dim>
void test()
{
    const int space_dim = dim + codim;
    using Function = functions::LinearFunction<dim, codim, space_dim>;
    using Mapping   = NewMapping<dim, codim>;

    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<space_dim; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    auto grid = CartesianGrid<dim>::create(3);
    auto F = Function::create(grid, A, b);

    auto flag = NewValueFlags::inv_hessian|NewValueFlags::inv_gradient|NewValueFlags::point;
    auto quad = QGauss<dim>(2);

    auto map = Mapping::create(F);
    map->template reset<sub_dim>(flag, quad);

    auto elem = map->begin();
    auto end  = map->end();
    const int s_id = 0;
    map->template init_cache<sub_dim>(elem);
    for (; elem != end; ++elem)
    {
        map->template fill_cache<sub_dim>(elem, s_id);
        elem->template get_inverse_values<1, sub_dim> (s_id).print_info(out);
        out << endl;
        elem->template get_inverse_values<2, sub_dim> (s_id).print_info(out);
        out << endl;
    }

}


int main()
{
    test<2,0>();
//    test<3,3>();

    return 0;
}

