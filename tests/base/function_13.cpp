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
 *  Test for SubFunction class
 *  author: pauletti
 *  date: Oct 12, 2014
 */

#include "../tests.h"

#include <igatools/functions/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/functions/sub_function.h>


template<int dim, int codim, int range>
void values_of_F(Function<dim, codim, range> &F)
{
    auto elem = F.begin();
    auto end  = F.end();

    const auto topology = Topology<dim>();

    F.init_cache(elem, topology);
    for (; elem != end; ++elem)
    {
        out.begin_item("Element " + std::to_string(elem->get_flat_index()));

        F.fill_cache(elem, topology, 0);
//        elem->get_points().print_info(out);
//        out << endl;
        out.begin_item("Value ");
        elem->template get_values<_Value, dim>(0).print_info(out);
        out.end_item();
//        elem->template get_values<1>().print_info(out);
//        out << endl;
//        elem->template get_values<2>().print_info(out);
//        out << endl;

        out.end_item();
    }
}


template<int dim, int codim, int range>
void create_fun()
{
    OUTSTART
    using Function = functions::LinearFunction<dim, codim, range>;

    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<range; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }

    auto flag = ValueFlags::point | ValueFlags::value |
                ValueFlags::gradient;

    auto grid = CartesianGrid<dim>::create(3);
    auto F = Function::create(grid, IdentityFunction<dim>::create(grid), A, b);


    const int k = dim - 1;
    for (auto &s_id : UnitElement<dim>::template elems_ids<k>())
    {
        using Grid =  CartesianGrid<dim>;
        using SubGridMap = typename Grid::template InterGridMap<k>;
        SubGridMap elem_map;
        auto sub_grid = grid->template get_sub_grid<k>(s_id, elem_map);

        auto subF = SubMapFunction<k, dim, range>::create(sub_grid, *F, s_id, elem_map);

        auto sub_quad = QGauss<k>(1);
        subF->reset(flag, sub_quad);
        out.begin_item("Face: " + std::to_string(s_id));
        values_of_F<k, 0, range>(*subF);
        out.end_item();
    }
    OUTEND
}


int main()
{
    create_fun<2, 0, 2>();
    return 0;
}

