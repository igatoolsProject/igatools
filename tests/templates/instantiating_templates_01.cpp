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
 * Test for instantiating templates
 * This test should fail if configure option was not set
 * author: pauletti
 * date: Jul 08, 2014
 *
 */

#include "../tests.h"
#include "igatools/base/function-template.h"


template <int dim, int space_dim>
class MyFun : public Function<dim, space_dim, 2>
{
public:
    using Base =  Function<dim, space_dim, 2>;
    using typename Base::Point;
    using typename Base::Value;
    using typename Base::Gradient;
    using typename Base::Hessian;


    void evaluate(const std::vector<Point> &points,
                  std::vector<Value> &values) const
    {}

    void evaluate_gradients(
        const std::vector<Point> &points,
        std::vector<Gradient> &gradient) const
    {}

    void evaluate_hessians(
        const std::vector<Point> &points,
        std::vector<Hessian> &hessians) const
    {}
};


int main()
{
    MyFun<2,3> f;
    return 0;
}
