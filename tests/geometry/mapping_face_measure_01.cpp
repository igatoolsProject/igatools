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
 *  Test for the spherical mapping class.
 *
 *  author: antolin
 *  date: 2014-03-19
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/function_element.h>
#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/mapping_element.h>

// TODO (pauletti, Nov 20, 2014): this test is more about cylindrical annulus
template <int dim, int sub_dim = dim-1>
void test_evaluate()
{
    using Function = functions::CylindricalAnnulus<dim>;
    using Mapping   = NewMapping<dim, 0>;

    auto grid = CartesianGrid<dim>::create();

    auto F = Function::create(grid, IdentityFunction<dim>::create(grid),
                              1., 2., 0., 1.0, 0.0, numbers::PI/2.0);
    Mapping map(F);

    auto flag = NewValueFlags::point | NewValueFlags::w_measure;
    auto quad = QGauss<sub_dim>(3);

    map.reset(flag, quad);

    auto elem = map.begin();
    auto end  = map.end();


    std::array<Real, UnitElement<dim>::template num_elem<sub_dim>()> face_area ;
    std::fill(face_area.begin(), face_area.end(), 0.0) ;

    map.template init_cache<sub_dim>(elem);

    for (; elem != end; ++elem)
    {
        if (elem->is_boundary())
        {
            for (auto &s_id : UnitElement<dim>::template elems_ids<sub_dim>())
            {
                if (elem->is_boundary(s_id))
                {
                    map.template fill_cache<sub_dim>(elem, s_id);
                    auto w_meas = elem->template get_w_measures<sub_dim>(s_id);
                    for (auto &w : w_meas)
                        face_area[s_id] += w;
                }
            }
        }
    }

    out << "Dimension " << dim << endl;
    for (Index face_id = 0; face_id < UnitElement<dim>::template num_elem<sub_dim>(); ++face_id)
    {
        out << "Area of face " << face_id << " : " << face_area[face_id] << endl;
    }


}

int main()
{
    out.depth_console(10);

//    test_evaluate<2>();
    test_evaluate<3>();

}
