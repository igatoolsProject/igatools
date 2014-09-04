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

#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/mapping_slice.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_element_accessor.h>

template <int dim,int codim>
void run_test()
{
    const int space_dim = dim + codim;
    auto grid = CartesianGrid<dim>::create();

    Derivatives<dim,space_dim,1, 1> A;
    Points<space_dim> b;
    //Dilation
    for (int i=0; i<dim; ++i)
        A[i][i] = i+1;
    //Traslation
    for (int i=0; i<dim; ++i)
        b[i] = i+1;

    auto map = LinearMapping<dim,codim>::create(grid, A, b);


    out << "Linear mapping" << "<" << dim << "," << space_dim << ">" << endl;
    out << "A =" << endl << A << endl;
    out << "b =" << b << endl << endl;

    QGauss<dim> quad(1);
    auto elem = map->begin();
    elem->init_cache(ValueFlags::point|ValueFlags::map_gradient, quad);

    elem->fill_cache();
    elem->get_map_values().print_info(out);
    elem->get_map_gradients().print_info(out);


    for (int face_id = 0; face_id < UnitElement<dim>::faces_per_element; ++face_id)
    {
        auto elem_map = std::make_shared<typename CartesianGrid<dim>::FaceGridMap>();
        auto face_grid = grid->get_face_grid(face_id, *elem_map);
        auto face_map =
            MappingSlice<dim-1,codim+1>::create(map, face_id, face_grid, elem_map);


        QGauss<dim-1> face_quad(1);
        auto face_elem = face_map->begin();
        face_elem->init_cache(ValueFlags::point|ValueFlags::map_gradient, face_quad);

        face_elem->fill_cache();

        out << "Map Values (x1,x2,...):" << endl;
        face_elem->get_map_values().print_info(out);
        out << endl;

        out << "Map Gradients (x1,x2,...):" << endl;
        face_elem->get_map_gradients().print_info(out);
        out << endl;
    }
}


int main()
{
    out.depth_console(10);

    run_test<2,0>();

    return  0;
}
