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

/**
 * This file contains the test for the Gauss-Legendre quadrature scheme,
 * implemented by the QGauss class.
 * \author Sebasian Pauletti (spauletti@gmail.com)
 * \author Massimiliano Martinelli (massimiliano.martinelli@gmail.com)
 * \date 08/04/2013
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>


void do_test_zero()
{
    out << "-------------------------------" << endl ;
    out << "Dimension: 0" << endl;
    QGauss<0> quad(1) ;
    quad.print_info(out) ;
    out << "-------------------------------" << endl ;
    out <<endl ;
}

//Isotropic
template <int dim>
void do_test()
{
    out << "-------------------------------" << endl ;
    const int max_deg = 4;
    out << "Dimension: " << dim << endl;
    for (int i = 1; i < max_deg; ++i)
    {
        out << "Points per direction: " << i << endl;
        QGauss<dim> quad(i);
        quad.print_info(out) ;
    }
    out << "-------------------------------" << endl ;
    out <<endl ;
}


//Anisotropic
template <int dim>
void do_test_aniso()
{
    out << "-------------------------------" << endl ;
    const int max_deg = 4;
    out << "Dimension: " << dim << endl;
    for (int i = 1; i < max_deg; ++i)
    {
        out << "Points per direction: ";
        TensorSize<dim> n_points;
        for (int j = 0; j < dim; ++j)
        {
            n_points[j] = i + j;
            out << i+j << " ";
        }
        out << endl;
        QGauss<dim> quad(n_points);
        quad.print_info(out) ;
    }
    out << "-------------------------------" << endl ;
    out <<endl ;
}


int main(int argc, char *argv[])
{
//    do_test_zero() ;

    do_test<1>();
    do_test<2>();
    do_test<3>();

    do_test_aniso<1>();
    do_test_aniso<2>();
    do_test_aniso<3>();

    return 0;
}
