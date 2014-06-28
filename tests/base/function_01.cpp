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
 *  Test for Function class, we define a linear function
 *  author: pauletti
 *  date: Jun 19, 2014
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_lib.h>

template<int dim>
class Func : public Function<dim>
{
    using typename Function<dim>::Point;
    using typename Function<dim>::Value;

public:
    void evaluate(const vector<Point> &points,
                  vector<Value> &values) const
    {
        const int NumPoints = points.size() ;
        for (int iPoint = 0 ; iPoint < NumPoints ; iPoint++)
        {
            values[ iPoint ] = 0.;
            for (int j = 0; j < dim; ++j)
                values[ iPoint ][0] += points[ iPoint ][j];
        }
    }
};



template <int dim>
void
run_test()
{
    const int n_pts = 3;
    out << "======  Testing functions f : R^" << dim << " --> R  ======" << endl;
    Func<dim> F ;
    QUniform<dim> quad(n_pts);
    auto p = quad.get_points().get_flat_cartesian_product();
    vector< typename Func<dim>::Value > v(p.size());

    F.evaluate(p,v);
    out << "Points = " << p << endl ;
    out << "function f(x) = sum x_i" << endl ;
    out << "Values= " << v << endl ;

    out << "===============================================" << endl ;

}



int main()
{
    run_test<0>();
    run_test<1>();
    run_test<2>();
    run_test<3>();

    return 0;
}
