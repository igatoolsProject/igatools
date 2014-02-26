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
 *  Testing quadratur extension from dim to dim+1
 *  author: antolin
 *  date: 2013-04-02
 *
 */

#include "../tests.h"
#include <igatools/base/quadrature.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/unit_element.h>
#include <igatools/base/types.h>




int main(int argc, char *argv[])
{


    const int num_pts = 3;

    out << "======  Extension from 2D to 3D  ======" << endl ;

    QGauss<2> quad_surf_2d(num_pts) ;

    out << "Points of the original quadrature" << endl ;
    out << quad_surf_2d.get_points() << endl ;
    out << "Weights of the original quadrature" << endl ;
    out << quad_surf_2d.get_weights() << endl ;
    out << endl ;


    for (int i = 0; i < UnitElement<3>::faces_per_element; ++i)
    {
        Quadrature<3> quad_surf_2d_3d = extend_quad_dim<2>(quad_surf_2d, i) ;
        out << endl ;
        out << "Points of the quadrature of the face " << i << endl ;
        out << quad_surf_2d_3d.get_points() << endl ;
        out << "Weights of the quadrature of the face " << i << endl ;
        out << quad_surf_2d_3d.get_weights() << endl ;
    }

    out << "=======================================" << endl ;
    out << endl ;
    out << endl ;

    out << "======  Extension from 1D to 2D  ======" << endl ;

    QGauss<1> quad_surf_1d(num_pts) ;

    out << "Points of the original quadrature" << endl ;
    out << quad_surf_1d.get_points() << endl ;
    out << "Weights of the original quadrature" << endl ;
    out << quad_surf_1d.get_weights() << endl ;
    out << endl ;


    for (int i = 0; i < UnitElement<2>::faces_per_element; ++i)
    {
        Quadrature<2> quad_surf_1d_2d = extend_quad_dim<1>(quad_surf_1d, i) ;
        out << endl ;
        out << "Points of the quadrature of the face " << i << endl ;
        out << quad_surf_1d_2d.get_points() << endl ;
        out << "Weights of the quadrature of the face " << i << endl ;
        out << quad_surf_1d_2d.get_weights() << endl ;
    }

    out << "=======================================" << endl ;


    return (0) ;
}

