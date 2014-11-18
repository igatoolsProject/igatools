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
 *  Test for the MappingSlice class
 *
 *  author: pauletti
 *  date: 2013-10-20
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_slice.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/function_element.h>

template <int sub_dim, int dim, int codim = 0>
void sub_function(const int n_knots = 3)
{
    using Func = IdentityFunction<dim>;
    using GridType = CartesianGrid<dim>;
    using SubFunc = SubFunction<sub_dim, dim, codim, dim, 1>;
    auto grid = GridType::create(n_knots);
    auto func = Func::create(grid);

    for (auto &s_id : UnitElement<dim>::template elems_ids<sub_dim>())
    {
        using  InterGridMap = typename GridType::template InterGridMap<sub_dim>;
        auto elem_map = std::make_shared<InterGridMap>(InterGridMap());

        auto sub_grid = grid->template get_sub_grid<sub_dim>(s_id, *elem_map);
        auto sub_func = SubFunc::create(sub_grid, func, s_id, *elem_map);


        out << "Face: " << s_id << endl;

        QGauss<sub_dim> quad(1);
        auto sub_func_flag = NewValueFlags::point | NewValueFlags::value;
        sub_func->reset(sub_func_flag, quad);

        auto f_elem =  sub_func->begin();
        auto end    =  sub_func->end();
        sub_func->init_cache(f_elem, Int<sub_dim>());


        for (; f_elem != end; ++f_elem)
        {
            out << "face element: " <<  f_elem->get_flat_index() << endl;
            sub_func->fill_cache(f_elem, 0, Int<sub_dim>());
            f_elem->template get_values<0,sub_dim>(0).print_info(out);
            out << endl;
        }

//        auto elem = face_map->begin();
//        auto end = face_map->end();
//        face_elem->init_cache(ValueFlags::point|ValueFlags::map_gradient, face_quad);
//        for (; face_elem != end; ++face_elem)
//        {
//            out << "face element: " <<  face_elem->get_flat_index() << endl;
//            face_elem->fill_cache();
//            face_elem->get_map_values().print_info(out);
//            out << endl;
//            //face_elem->get_gradients().print_info(out);
//        }
        out << endl;
    }
}


int main()
{
    out.depth_console(10);

    sub_function<1, 2>();

    return  0;
}
