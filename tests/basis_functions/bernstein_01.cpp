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
 *  @brief Test for bernstein polynomial class
 *  @author pauletti
 *  @date 2013-01-10
 */

#include "../tests.h"

#include <igatools/basis_functions/bernstein_basis.h>
#include <boost/numeric/ublas/io.hpp>

int main()
{
    out.depth_console(10);

    const int n=5;
    SafeSTLVector< Real > points(n);
    for (int i=0; i<n; ++i)
        points[i] = Real(i)/(n-1);
    out << "points: ";
    points.print_info(out);
    out << endl;

    for (int p = 0; p<3; p++)
    {
        out << "degree: " << p << endl;

        auto values = BernsteinBasis::derivative(0, p, points);
        out << "values: " << values << endl;

        values = BernsteinBasis::derivative(1, p, points);
        out << "first derivatives: " << values << endl;


        values = BernsteinBasis::derivative(2, p, points);
        out << " second derivatives: " << values << endl;

        out <<  endl;
    }

    return 0;
}
