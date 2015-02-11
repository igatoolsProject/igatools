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
 *  Test for the  CartesianGrid::get_face_grid
 *
 *  author: pauletti
 *  date: 2013-01-27
 *  QA: v0.2 (pauletti, 2013-10-25)
 */

#include "../tests.h"
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>

template<int dim>
void face_uniform(const int n_knots)
{
    auto grid = CartesianGrid<dim>::create(n_knots);
    out.begin_item("CaresianGrid<" + std::to_string(dim) +">");
    grid->print_info(out);
    out.end_item();
    for (int i = 0; i < UnitElement<dim>::faces_per_element; ++i)
    {
        out << "Face: " << i << endl;
        std::map<int, int> elem_map;
        auto face_grid = grid->get_face_grid(i, elem_map);
        face_grid->print_info(out);

        out << "Face-grid to element-grid map:" << endl;
        for (auto x : elem_map)
            out << x.first << " " << x.second << endl;
        out << endl;
    }
    out << endl;
}


template<int dim>
void face_non_uniform()
{
    using Grid =  CartesianGrid<dim>;
    TensorSize<dim> n_knots;
    for (int i = 0; i < dim; ++i)
        n_knots[i] = 2+i;
    auto grid = Grid::create(n_knots);
    out.begin_item("CaresianGrid<" + std::to_string(dim) +">");
    grid->print_info(out);
    out.end_item();
    for (int i = 0; i < UnitElement<dim>::faces_per_element; ++i)
    {
        out << "Face: " << i << endl;
        typename Grid::FaceGridMap elem_map;
        auto face_grid = grid->get_face_grid(i, elem_map);
        face_grid->print_info(out);

        out << "Face-grid to element-grid map:" << endl;
        for (auto x : elem_map)
            out << x.first->get_flat_index() << " " << x.second->get_flat_index() << endl;
        out << endl;
    }
    out << endl;
}



int main()
{
    out.depth_console(10);

    face_non_uniform<1>();
    face_non_uniform<2>();
    face_non_uniform<3>();

    return  0;
}
