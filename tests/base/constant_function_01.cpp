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
 *  Test for the ConstantFunction.
 *
 *  author: pauletti
 *  date: 2013-11-01
 *
 */

#include "../tests.h"
#include <igatools/base/function_lib.h>

using namespace functions;

template<int dim, int rdim, int rank>
void
run_test(ConstantFunction<dim,rdim,rank> &fun)
{
    using value_t = typename Function<dim,rdim,rank>::Value;
    using point_t = typename Function<dim,rdim,rank>::Point;


    const int n_pts = 3;
    vector<point_t> points(n_pts);
    vector<value_t> values(n_pts);

    fun.evaluate(points, values);
    out << points << endl;
    out << values << endl;
}



int main()
{
    ConstantFunction<0> f0({Real(2.)});
    run_test<0,1,1>(f0);

    ConstantFunction<1> f1({Real(2.)});
    run_test<1,1,1>(f1);

    ConstantFunction<2> f2({Real(2.)});
    run_test<2,1,1>(f2);

    ConstantFunction<2,2> f22({Real(2.), Real(3.)});
    run_test<2,2,1>(f22);

    return 0;
}

