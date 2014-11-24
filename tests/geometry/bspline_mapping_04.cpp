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
 *  Test for IgMapping class.
 *
 *  author: martinelli
 *  date: 18/04/2013
 *
 */

// TODO (pauletti, Nov 20, 2014): the relation with the other test should be
// clarify

#include "../tests.h"

#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/mapping_element.h>
#include <igatools/base/ig_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/new_bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/function_element.h>



template <int dim, int codim = 0, int rank = 1>
void bspline_map(const int deg = 2)
{
    const int sub_dim = dim;

    using Space = NewBSplineSpace<dim, dim+codim, rank>;
    using Function = IgFunction<Space>;
    using Mapping   = NewMapping<dim, codim>;



    //----------------------------------------------------------------------------------------------
    out << "Dim: " << dim << endl ;
    int n_knots = 2;
    CartesianProductArray<Real , dim> coord ;
    for (int i = 0; i < dim; ++i)
    {
        vector<Real> tmp_coord;
        for (int j = 0; j < n_knots; ++j)
            tmp_coord.push_back(j);
        coord.copy_data_direction(i,tmp_coord);
    }


    auto grid = CartesianGrid<dim>::create(coord);
    auto space = Space::create(deg, grid);

    typename Function::CoeffType control_pts(space->get_num_basis());

    if (dim == 2)
    {
        int id = 0 ;

        // x coords
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 1.5 ;
        control_pts(id++) = 1.5 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 2.0 ;
        control_pts(id++) = 2.0 ;


        // y coords
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 1.5 ;
        control_pts(id++) = 1.5 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 2.0 ;
        control_pts(id++) = 2.0 ;
        control_pts(id++) = 0.0 ;
    }
    else if (dim == 3)
    {
        int id = 0 ;

        // x coords
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 1.5 ;
        control_pts(id++) = 1.5 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 2.0 ;
        control_pts(id++) = 2.0 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 1.5 ;
        control_pts(id++) = 1.5 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 2.0 ;
        control_pts(id++) = 2.0 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 1.5 ;
        control_pts(id++) = 1.5 ;

        control_pts(id++) = 0.0 ;
        control_pts(id++) = 2.0 ;
        control_pts(id++) = 2.0 ;


        // y coords
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 1.5 ;
        control_pts(id++) = 1.5 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 2.0 ;
        control_pts(id++) = 2.0 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 1.5 ;
        control_pts(id++) = 1.5 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 2.0 ;
        control_pts(id++) = 2.0 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 1.5 ;
        control_pts(id++) = 1.5 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 2.0 ;
        control_pts(id++) = 2.0 ;
        control_pts(id++) = 0.0 ;

        // z coords
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 0.0 ;
        control_pts(id++) = 0.0 ;

        control_pts(id++) = 0.5 ;
        control_pts(id++) = 0.5 ;
        control_pts(id++) = 0.5 ;
        control_pts(id++) = 0.5 ;
        control_pts(id++) = 0.5 ;
        control_pts(id++) = 0.5 ;
        control_pts(id++) = 0.5 ;
        control_pts(id++) = 0.5 ;
        control_pts(id++) = 0.5 ;

        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;
        control_pts(id++) = 1.0 ;

    }


    auto F = Function::create(space, control_pts);
    auto map = Mapping::create(F);

    auto quad = QGauss<dim>(3);
    auto flag =  NewValueFlags::value| NewValueFlags::gradient
                 | NewValueFlags::hessian;

    map->template reset<sub_dim>(flag, quad);

    auto elem = map->begin();
    auto end  = map->end();
    const int s_id = 0;

    map->template init_cache<sub_dim>(elem);
    for (; elem != end; ++elem)
    {
        map->template fill_cache<sub_dim>(elem, s_id);
        out << "Values : ";
        elem->template get_values<0,sub_dim>(s_id).print_info(out);
        out << endl;
        out << "Gradients : ";
        elem->template get_values<1,sub_dim>(s_id).print_info(out);
        out << endl;
        out << "Hessians : ";
        elem->template get_values<2,sub_dim>(s_id).print_info(out);
        out << endl;
    }


}

int main()
{
    out.depth_console(10);

    bspline_map<2>();
    bspline_map<3>();

}
