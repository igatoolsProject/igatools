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

/**
 *  @file
 *  @brief  get_weihts()
 *  @author pauletti
 *  @date   2015-08-19
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/grid_cache_handler.h>
#include <igatools/geometry/grid_element.h>


template <int dim, int sdim = dim>
void elem_weights(const int n_knots = 5)
{
    OUTSTART

    using Grid = CartesianGrid<dim>;
    using Flags = typename Grid::ElementAccessor::Flags;
    auto grid = Grid::create(n_knots);

    auto flag = Flags::w_measure;
    auto cache_handler = grid->create_cache_handler();
    cache_handler->template set_flags<sdim>(flag);

    auto quad = QGauss<sdim>::create(2);
    auto elem = grid->begin();
    cache_handler->template init_cache<sdim>(elem, quad);

    for (; elem != grid->end(); ++elem)
    {
        for (auto &s_id : UnitElement<dim>::template elems_ids<sdim>())
        {
            cache_handler->template fill_cache<sdim>(elem, s_id);
            elem->template get_w_measures<sdim>(s_id).print_info(out);
            out << endl;
        }
        out << endl;
    }

    OUTEND
}



int main()
{
    elem_weights<0>();
    elem_weights<1>();
    elem_weights<2>();
    elem_weights<3>();
    elem_weights<1,0>();
    elem_weights<2,1>();
    elem_weights<3,2>();

    return  0;
}
