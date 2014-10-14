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
 *  Development of AnaliticalFunction
 *  author: pauletti
 *  date: Oct 12, 2014
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/formula_function.h>
#include <igatools/base/function_element.h>


template<int dim, int range>
class LinearFunction : public FormulaFunction<dim, 0, range, 1>
{
public:
    using parent_t = FormulaFunction<dim, 0, range, 1>;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    LinearFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
                   const ValueFlags &flag, const Quadrature<dim> &quad,
                   const Gradient &A, const Value &b)
        :
        parent_t::FormulaFunction(grid, flag, quad),
        A_(A),
        b_(b)
    {}


private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const
    {
        auto point = points.begin();
        for (auto &val : values)
        {
            val = action(A_, *point) + b_;
            ++point;
        }
    }

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const
    {
        for (auto &val : values)
            val = A_;
    }

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const
    {
        for (auto &val : values)
            val = 0.;
    }

private:
    const Gradient A_;
    const Value    b_;
};



template<int dim, int range>
void test(NewFunction<dim, 0, range> &F, shared_ptr<CartesianGrid<dim>> grid)
{
    GridForwardIterator<FunctionElement<dim, 0, range,1>> elem(grid, 0);
    GridForwardIterator<FunctionElement<dim, 0, range,1>> end(grid, IteratorState::pass_the_end);

    F.init_element(elem);
    for (; elem != end; ++elem)
    {
        F.fill_element(elem);
        elem->get_points().print_info(out);
        out << endl;
        elem->get_values().print_info(out);
        out << endl;
        elem->get_gradients().print_info(out);
        out << endl;
        elem->get_hessians().print_info(out);
        out << endl;
    }
}


template<int dim, int range>
void create_fun()
{
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

    auto flag = ValueFlags::point | ValueFlags::value | ValueFlags::gradient |
                ValueFlags::hessian;
    auto quad = QGauss<dim>(2);
    auto grid = CartesianGrid<dim>::create(3);
    Function F(grid, flag, quad, A, b);
    test<dim,range>(F, grid);
}




int main()
{
    create_fun<2,2>();
//    test<3,3>();

    return 0;
}

