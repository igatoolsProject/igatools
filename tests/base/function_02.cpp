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
 *  Test for user defined linear function based on the virtual class Function.
 */

#include "../tests.h"

#include <igatools/base/function.h>

template<int dim, int range>
class LinearFunction : public Function<dim, range, 1 >
{


public:
    LinearFunction()
    {

        for (int i=0; i<range; i++)
        {
            for (int j=0; j<dim; j++)
                if (j == i)
                    A[j][j] = 2.;

            b[i] = i;
        }

    }

    void evaluate(
        const std::vector < typename LinearFunction<dim, range>::PointType >   &points,
              std::vector < typename LinearFunction<dim, range>::ValueType >   &values) const
    {

        for (int i=0; i<points.size(); i++)
            values[i] = action(A,points[i]) + b;
    }

private:
    typename LinearFunction<dim, range>::GradientType A;
    typename LinearFunction<dim, range>::ValueType    b;

};



int main()
{
    const int dim=2;
    const int range=2;

    LinearFunction< dim, range > F;

    std::vector < LinearFunction< dim, range >::PointType > x(2);
    std::vector < LinearFunction< dim, range >::ValueType > y(2);
    x[1][0] = 1.;

    F.evaluate(x,y);

    out << x << endl;
    out << y << endl;

    return 0;
}

