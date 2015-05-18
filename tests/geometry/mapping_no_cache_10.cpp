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
 *  Test for IgFunction in a Mapping without the use of the cache
 *  author: martinelli
 *  date: Jan 30, 2015
 */

#include "../tests.h"

#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>
#include <igatools/base/ig_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/function_element.h>


template<int dim>
void ig_mapping(const int deg = 1)
{
    OUTSTART

    using Space = BSplineSpace<dim,dim>;
    using RefSpace = ReferenceSpace<dim, dim>;
    using Function = IgFunction<dim,0,dim,1>;
    using Mapping   = Mapping<dim,0>;


    auto quad = QGauss<dim>(2);
    auto grid = CartesianGrid<dim>::create(3);

    auto space = Space::create(deg, grid);
    auto coeff = EpetraTools::create_vector(space, DofProperties::active);
    (*coeff)[0] = 1.;
    auto F = Function::create(space, coeff);

    auto map = Mapping::create(F);

    auto elem = map->begin();
    auto end  = map->end();

    for (; elem != end; ++elem)
    {
        elem->template evaluate_at_points<_Value>(quad).print_info(out);
        out << endl;
        elem->template evaluate_at_points<_Gradient>(quad).print_info(out);
        out << endl;
        elem->template evaluate_at_points<_Hessian>(quad).print_info(out);
        out << endl;
    }

    OUTEND
}


int main()
{
    ig_mapping<2>();
    ig_mapping<3>();

    return 0;
}

