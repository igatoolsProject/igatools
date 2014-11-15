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
 *  Test for the l2_projection function.
 *  Bspline spaces case
 *
 *  author: pauletti
 *  date: 2013-10-10
 *  QA: The 0 dim case not returning appropriate value
 */


#include <igatools/base/quadrature_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/base/formula_function.h>



// TODO (pauletti, Nov 13, 2014): delete this after cleaning and arranging test
using numbers::PI;

template<int dim>
class BoundaryFunction : public FormulaFunction<dim>
{
private:
   using base_t = NewFunction<dim>;
    using parent_t = FormulaFunction<dim>;
   using self_t = BoundaryFunction<dim>;
    using typename base_t::GridType;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    template <int order>
        using Derivative = typename parent_t::template Derivative<order>;
public:
    BoundaryFunction(std::shared_ptr<GridType> grid)
    : FormulaFunction<dim>(grid, IdentityFunction<dim>::create(grid))
      {}

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid)
    {
        return std::shared_ptr<base_t>(new self_t(grid));
    }

    std::shared_ptr<base_t> clone() const override
    {
        return std::make_shared<self_t>(self_t(*this));
    }

    Real value(Points<dim> x) const
    {
        Real f = 1;
        for (int i = 0; i<dim; ++i)
            f = f * cos(2*PI*x[i]);
        return f;
    }

    void evaluate_0(const ValueVector<Point> &points,
                       ValueVector<Value> &values) const override
    {
        for (int i = 0; i<points.size(); ++i)
        {
            Points<dim> p = points[i];
            values[i][0] = this->value(p);
        }
    }
    void evaluate_1(const ValueVector<Point> &points,
                       ValueVector<Derivative<1>> &values) const override
                               {}

       void evaluate_2(const ValueVector<Point> &points,
                       ValueVector<Derivative<2>> &values) const override
                               {}
};




