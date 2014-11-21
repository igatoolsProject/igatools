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
 *  Test of for the faces normals using a cylindrical annulus mapping.
 *  author: antolin
 *  date: 2014-03-18
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/function_element.h>
#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/mapping_element.h>


template <int dim>
auto create_mapping1(shared_ptr<const CartesianGrid<dim>> grid)
{
    using Function = functions::CylindricalAnnulus<dim>;

    return
        Function::create(grid, IdentityFunction<dim>::create(grid),
                         1, 2, 0, 2.0, 0.0, numbers::PI/2.0);
}

template <int dim>
auto create_mapping2(shared_ptr<const CartesianGrid<dim>> grid)
{
    using Function = IdentityFunction<dim>;

    return Function::create(grid);
}

template <int sub_dim, int dim, int codim =0 >
void boundary_normals()
{
    using Mapping = NewMapping<dim, codim>;

    auto grid = CartesianGrid<dim>::create();
    auto map_func =  create_mapping1<dim>(grid);

    auto flag = NewValueFlags::w_measure|NewValueFlags::point|NewValueFlags::outer_normal;
    auto quad = QGauss<sub_dim>(1);

    Mapping map(map_func);
    map.reset(flag, quad);

    auto elem = map.begin();
    auto end = map.end();

    map.template init_cache<sub_dim>(elem);
    for (; elem != end; ++elem)
    {
        for (auto &s_id : UnitElement<dim>::template elems_ids<sub_dim>())
        {
            out << "Face: " << s_id << endl;
            out << "  Normal vector:" << endl;
            map.template fill_cache<sub_dim>(elem, s_id);
            elem->template get_boundary_normals<sub_dim>(s_id).print_info(out);
            out << endl;
        }
    }



}

int main()
{
    boundary_normals<2,3,0>();
    return 0;
}

