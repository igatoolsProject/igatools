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
 *  Test for MappingSlice class
 *
 *  author: pauletti
 *  date: Feb 2, 2013
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/sub_function.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/function_element.h>
#include <igatools/base/identity_function.h>


template <int sub_dim, int dim, int space_dim>
void sub_map(const int n_knots = 2)
{
    OUTSTART
    using Func = functions::LinearFunction<dim, 0, space_dim>;
    using GridType = CartesianGrid<dim>;
    using SubFunc = SubMapFunction<sub_dim, dim, space_dim>;

    auto grid = GridType::create(n_knots);

    typename Func::Gradient A;
    typename Func::Value    b;
    for (int i=0; i<dim; ++i)
        A[i][i] = i+1;
    for (int i=0; i<dim; ++i)
        b[i] = i+1;

    auto func = Func::create(grid, IdentityFunction<dim>::create(grid), A, b);

    for (auto &s_id : UnitElement<dim>::template elems_ids<sub_dim>())
    {
        using  InterGridMap = typename GridType::template InterGridMap<sub_dim>;
        auto elem_map = std::make_shared<InterGridMap>(InterGridMap());

        auto sub_grid = grid->template get_sub_grid<sub_dim>(s_id, *elem_map);
        auto sub_func = SubFunc::create(sub_grid, func, s_id, *elem_map);

        QGauss<sub_dim> f_quad(1);
        auto sub_func_flag = ValueFlags::point | ValueFlags::value
                             | ValueFlags::gradient;
        sub_func->reset(sub_func_flag, f_quad);

        auto f_elem =  sub_func->begin();
        auto end    =  sub_func->end();
        sub_func->init_cache(f_elem, Int<sub_dim>());


        for (; f_elem != end; ++f_elem)
        {
            sub_func->fill_cache(f_elem, 0, Int<sub_dim>());
            f_elem->template get_values<0,sub_dim>(0).print_info(out);
            out << endl;
            f_elem->template get_values<1,sub_dim>(0).print_info(out);
            out << endl;
        }
    }
    OUTEND
}


int main()
{
    out.depth_console(10);

    sub_map<1, 2, 2>();

    return  0;
}
