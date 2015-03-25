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
 *  Test for IgFunction in a Mapping
 *  author: pauletti
 *  date: Oct 23, 2014
 */

#include "../tests.h"

#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>
#include <igatools/base/ig_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/function_element.h>


template<int dim, int codim, int sub_dim=dim>
void ig_mapping(const int deg = 1)
{
    OUTSTART

    using Space = BSplineSpace<dim, dim+codim>;
    using RefSpace = ReferenceSpace<dim, dim+codim>;
    using Function = IgFunction<RefSpace>;
    using Mapping   = Mapping<dim, codim>;


    auto flag =  ValueFlags::value| ValueFlags::gradient
                 | ValueFlags::hessian;
    auto quad = QGauss<dim>(2);
    auto grid = CartesianGrid<dim>::create(3);

    auto space = Space::create(deg, grid);
    vector<Real> coeff(space->get_num_basis());
    coeff[0] = 1.;
    auto F = Function::create(space, IgCoefficients(*space,DofProperties::active,coeff));

    auto map = Mapping::create(F);
    map->template reset<sub_dim>(flag, quad);

    auto elem = map->begin();
    auto end  = map->end();
    const int s_id = 0;

    map->template init_cache<sub_dim>(elem);
    for (; elem != end; ++elem)
    {
        map->template fill_cache<sub_dim>(elem, s_id);

        elem->template get_values<0,sub_dim>(s_id).print_info(out);
        out << endl;
        elem->template get_values<1,sub_dim>(s_id).print_info(out);
        out << endl;
        elem->template get_values<2,sub_dim>(s_id).print_info(out);
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

