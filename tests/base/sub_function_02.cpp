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
 *  Test for MappingSlice class
 *
 *  author: pauletti
 *  date: Feb 2, 2013
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_slice.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/function_element.h>
#include <igatools/base/identity_function.h>


template <int sub_dim, int dim, int space_dim>
void sub_function(const int n_knots = 2)
{
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

    out << "Linear mapping" << "<" << dim << "," << space_dim << ">" << endl;
    out << "A =" << endl << A << endl;
    out << "b =" << b << endl << endl;



    QGauss<dim> quad(1);
    auto func_flag = NewValueFlags::point | NewValueFlags::value
            | NewValueFlags::gradient;
    func->reset(func_flag, quad);
    auto elem =  func->begin();
    func->init_cache(elem, Int<dim>());
    func->fill_cache(elem, 0, Int<dim>());
    elem->template get_values<0,dim>(0).print_info(out);
    elem->template get_values<1,dim>(0).print_info(out);



    for (auto &s_id : UnitElement<dim>::template elems_ids<sub_dim>())
    {
        using  InterGridMap = typename GridType::template InterGridMap<sub_dim>;
        auto elem_map = std::make_shared<InterGridMap>(InterGridMap());

        auto sub_grid = grid->template get_sub_grid<sub_dim>(s_id, *elem_map);
        auto sub_func = SubFunc::create(sub_grid, func, s_id, *elem_map);


        //out << "Face: " << s_id << endl;

        QGauss<sub_dim> f_quad(1);
        auto sub_func_flag = NewValueFlags::point | NewValueFlags::value
                | NewValueFlags::gradient;
        sub_func->reset(sub_func_flag, f_quad);

        auto f_elem =  sub_func->begin();
        auto end    =  sub_func->end();
        sub_func->init_cache(f_elem, Int<sub_dim>());


        for (; f_elem != end; ++f_elem)
        {
           // out << "face element: " <<  f_elem->get_flat_index() << endl;
            sub_func->fill_cache(f_elem, 0, Int<sub_dim>());
            out << "Map Values (x1,x2,...):" << endl;
            f_elem->template get_values<0,sub_dim>(0).print_info(out);
            out << endl;
            out << "Map Gradients (x1,x2,...):" << endl;
            f_elem->template get_values<1,sub_dim>(0).print_info(out);
            out << endl;
        }
    }

//    for (int face_id = 0; face_id < UnitElement<dim>::n_faces; ++face_id)
//    {
//        auto elem_map = std::make_shared<typename CartesianGrid<dim>::FaceGridMap>();
//        auto face_grid = grid->get_face_grid(face_id, *elem_map);
//        auto face_map =
//            MappingSlice<dim-1,codim+1>::create(map, face_id, face_grid, elem_map);
//
//
//        QGauss<dim-1> face_quad(1);
//        auto face_elem = face_map->begin();
//        face_elem->init_cache(ValueFlags::point|ValueFlags::map_gradient, face_quad);
//
//        face_elem->fill_cache();
//
//        out << "Map Values (x1,x2,...):" << endl;
//        face_elem->get_map_values().print_info(out);
//        out << endl;
//
//        out << "Map Gradients (x1,x2,...):" << endl;
//        face_elem->get_map_gradients().print_info(out);
//        out << endl;
//    }
}


int main()
{
    out.depth_console(10);

    sub_function<1, 2, 2>();

    return  0;
}
