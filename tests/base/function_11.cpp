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
 *  Test for Function class, as a prototype for an analytical function
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_element.h>

template<int dim, int range>
class LinearFunction : public NewFunction<dim, 0, range>
{
public:
    using parent_t = NewFunction<dim, 0, range>;
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
        parent_t::NewFunction(grid, flag, quad),
        flag_(flag),
        quad_(quad),
        A_(A),
        b_(b)
    {}

    void init_element(ElementIterator &elem)
    {
        auto &el    = elem.get_accessor();
        GridUniformQuadCache<dim>::init_element_cache(el);
        auto &cache = this->get_cache(elem);
        if (cache == nullptr)
        {
            using Cache = typename ElementAccessor::CacheType;
            cache = shared_ptr<Cache>(new Cache);
        }
        cache->resize(flag_, quad_.get_num_points());
    }

    void fill_element(ElementIterator &elem)
    {
        auto &el    = elem.get_accessor();
        GridUniformQuadCache<dim>::fill_element_cache(el);
        const auto points = el.get_points();
        auto &cache = this->get_cache(elem);
        if (flag_.fill_values())
            evaluate_0(points, cache->values_);
        if (flag_.fill_gradients())
            evaluate_1(points, std::get<1>(cache->derivatives_));
        if (flag_.fill_hessians())
            evaluate_2(points, std::get<2>(cache->derivatives_));

    }

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
    ValueFlagsHandler flag_;
    Quadrature<dim> quad_;
    const Gradient A_;
    const Value    b_;
};




template<int dim, int range>
void test()
{
    //const int n_points = 2;
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


int main()
{
    test<2,2>();
//    test<3,3>();

    return 0;
}

