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
 *  Test for the CartesianGrid element iterator
 *  when getting face related values.
 *
 *  author: pauletti
 *  date: Aug 28, 2014
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/geometry/cartesian_grid_element.h>

template <int dim, int k = dim-1>
void face_values(const TensorSize<dim> &n_knots)
{
    OUTSTART
    using Grid = CartesianGrid<dim>;
    using ElementHandler = typename Grid::ElementHandler;

    auto grid = Grid::create(n_knots);

    NewValueFlags flag = NewValueFlags::measure|
                         NewValueFlags::w_measure|
                         NewValueFlags::point;
    QUniform<k> quad(2);
    ElementHandler cache(grid);
    cache.template reset<k>(flag, quad);

    auto elem = grid->begin();
    cache.template init_cache<k>(elem);

    for (; elem != grid->end(); ++elem)
    {
        out << "Element: ";
        elem->print_info(out);
        out << endl;

        for (auto &face_id : UnitElement<dim>::template elems_ids<k>())
        {
            if (elem->is_boundary(face_id))
            {
                cache.template fill_cache<k>(elem, face_id);
                out << "face: " << face_id << endl;

                out << "meas: "<< elem->template get_measure<k>(face_id) << endl;

                out << "w_meas: " << endl;
                elem->template get_w_measures<k>(face_id).print_info(out);
                out << endl;

                out << "points: " << endl;
                elem->template get_points<k>(face_id).print_info(out);
                out << endl;
            }
        }
    }

    OUTEND
}


int main()
{
    out.depth_console(10);

    face_values(TensorSize<1>(3));
    face_values(TensorSize<2>(3));
    face_values(TensorSize<3>(3));

    face_values(TensorSize<1>(arr::sequence<1>(2)));
    face_values(TensorSize<2>(arr::sequence<2>(2)));
    face_values(TensorSize<3>(arr::sequence<3>(2)));

    return  0;
}
