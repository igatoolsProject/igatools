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

#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/mapping_slice.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_element_accessor.h>

template <int dim,int codim>
void run_test()
{
    //const int space_dim = dim + codim;
    const int n_elem = 2;
    auto grid = CartesianGrid<dim>::create(n_elem+1);
    auto map = IdentityMapping<dim>::create(grid);

    for (int face_id = 0; face_id < UnitElement<dim>::faces_per_element; ++face_id)
    {
        out << "Face: " << face_id << endl;
        auto elem_map = std::make_shared<typename CartesianGrid<dim>::FaceGridMap>();
        auto face_grid = grid->get_face_grid(face_id, *elem_map);
        auto face_map =
            MappingSlice<dim-1,codim+1>::create(map, face_id, face_grid, elem_map);

        QGauss<dim-1> face_quad(1);
        auto face_elem = face_map->begin();
        auto end = face_map->end();
        face_elem->init_cache(ValueFlags::point|ValueFlags::map_gradient, face_quad);
        for (; face_elem != end; ++face_elem)
        {
            out << "face element: " <<  face_elem->get_flat_index() << endl;
            face_elem->fill_cache();
            face_elem->get_map_values().print_info(out);
            //face_elem->get_gradients().print_info(out);
        }
        out << endl;
    }
}


int main()
{
    out.depth_console(10);

    run_test<2,0>();

    return  0;
}
