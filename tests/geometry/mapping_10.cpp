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
 *  Test for IgFunction in a Mapping
 *  author: pauletti
 *  date: Oct 23, 2014
 */

#include "../tests.h"

#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/mapping_element.h>
#include <igatools/../../source/geometry/grid_forward_iterator.cpp>
#include <igatools/base/ig_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/new_bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/function_element.h>


template<int dim, int codim>
void ig_mapping(const int deg = 1)
{
    OUTSTART

    using Space = NewBSplineSpace<dim, dim+codim>;
    using Function = IgFunction<Space>;
    using Mapping   = NewMapping<dim, codim>;
    using ElementIt = typename Mapping::ElementIterator;

    auto flag =  NewValueFlags::value| NewValueFlags::gradient
            | NewValueFlags::hessian;
    auto quad = QGauss<dim>(2);
    auto grid = CartesianGrid<dim>::create(3);

    auto space = Space::create(deg, grid);
    typename Function::CoeffType coeff(space->get_num_basis());
    coeff(0) = 1.;
    auto F = Function::create(flag, quad, space, coeff);

    Mapping map(F, flag, quad);
    ElementIt elem(grid, 0);
    ElementIt end(grid, IteratorState::pass_the_end);

    map.init_element(elem);
    for (; elem != end; ++elem)
    {
        map.fill_element(elem);
        elem->get_values().print_info(out);
        out << endl;
        elem->get_gradients().print_info(out);
        out << endl;
        elem->get_hessians().print_info(out);
        out << endl;
    }

    OUTEND
}


int main()
{
    ig_mapping<2,0>();
    ig_mapping<3,0>();

    return 0;
}

