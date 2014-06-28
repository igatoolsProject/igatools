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
#include <igatools/base/function.h>

template<int dim, int range>
class LinearFunction : public Function<dim, range>
{
public:
    using typename Function<dim, range>::Point;
    using typename Function<dim, range>::Value;
    using typename Function<dim, range>::Gradient;

    LinearFunction(const Gradient &A, const Value &b)
        :
        A_ {A},
       b_ {b}
    {}

    void evaluate(const vector<Point> &points,
                  vector<Value> &values) const
    {
        auto point = points.begin();
        for (auto &val : values)
        {
            val = action(A_, *point) + b_;
            ++point;
        }
    }

private:
    const Gradient A_;
    const Value    b_;
};



template<int dim, int range>
void test()
{
    const int n_points = 2;
    using Function = LinearFunction<dim, range>;

    typename Function::Value    b;
    typename Function::Gradient A;
    for (int i=0; i<range; i++)
    {
        for (int j=0; j<dim; j++)
            if (j == i)
                A[j][j] = 2.;
        b[i] = i;
    }


    Function F(A,b);
    vector<typename Function::Point> x(n_points);
    vector<typename Function::Value> y(n_points);
    x[1][0] = 1.;

    F.evaluate(x,y);

    out << x << endl;
    out << y << endl;
}


int main()
{
    test<2,2>();
    test<3,3>();

    return 0;
}

