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
#include <igatools/base/new_function.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/../../source/geometry/grid_forward_iterator.cpp>
IGA_NAMESPACE_OPEN
template<int dim, int range = 1, int rank = 1>
class FunctionElement : public CartesianGridElement<dim>
{
public:
    using Func = NewFunction<dim, range, rank>;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Gradient = typename Func::Gradient;
    using Hessian  = typename Func::Hessian;
    using ContainerType = CartesianGrid<dim>;

private:
//  template<int k>
//  ValueVector<Derivative<k>> get_derivative() const;
    template <int order>
    using Derivative = typename Func::template Derivative<order>;
public:
    using CartesianGridElement<dim>::CartesianGridElement;
    ValueVector<Point> get_points() const
    {
        return CartesianGridElement<dim>::get_points();
    }

    ValueVector<Value> const &get_values() const
    {
        return elem_cache_->values_;
    }

    template<int order>
    auto const &get_derivative() const
    {
        return std::get<order>(elem_cache_->derivatives_);
    }

    ValueVector<Gradient> const &get_gradients() const
    {
        return get_derivative<1>();
    }

    ValueVector<Hessian> const &get_hessians() const
    {
        return get_derivative<2>();
    }

private:
    struct Cache : public CacheStatus
    {
        void resize(const ValueFlagsHandler &flags_handler,
                    const int n_points)
        {
            //TODO(pauletti, Oct 11, 2014): missing all necesary clears
            flags_handler_ = flags_handler;

            if (flags_handler_.fill_values())
                values_.resize(n_points);

            if (flags_handler_.fill_gradients())
                std::get<1>(derivatives_).resize(n_points);

            if (flags_handler_.fill_hessians())
                std::get<2>(derivatives_).resize(n_points);

            set_initialized(true);
        }

        void print_info(LogStream &out) const
        {
            flags_handler_.print_info(out);
            values_.print_info(out);
            std::get<1>(derivatives_).print_info(out);
            std::get<2>(derivatives_).print_info(out);
        }

        ValueVector<Value> values_;
        std::tuple<ValueVector<Derivative<0>>,
            ValueVector<Derivative<1>>,
            ValueVector<Derivative<2>>> derivatives_;
        ValueFlagsHandler flags_handler_;
    };

    std::shared_ptr<Cache> elem_cache_;
public:
    using CacheType = Cache;
private:
    template <typename Accessor> friend class GridForwardIterator;
    friend class NewFunction<dim, range, rank>;
};

template<int dim, int range, int rank>
auto
NewFunction<dim,range,rank>::
get_cache(NewFunction<dim,range,rank>::ElementIterator &elem) -> std::shared_ptr<typename ElementAccessor::CacheType> &
{
    return elem.get_accessor().elem_cache_;
}
template class NewFunction<2,2,1>;
template class FunctionElement<2,2,1>;
template class GridForwardIterator<FunctionElement<2,2,1>>;
IGA_NAMESPACE_CLOSE

//template class CartesianGridElement<1>;
//template class GridForwardIterator<CartesianGridElement<1>>;

template<int dim, int range>
class LinearFunction : public NewFunction<dim, range>
{
public:
    using parent_t = NewFunction<dim, range>;
    using typename NewFunction<dim, range>::Point;
    using typename NewFunction<dim, range>::Value;
    using typename NewFunction<dim, range>::Gradient;
    using typename NewFunction<dim, range>::ElementIterator;
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
        A_ {A},
       b_ {b}
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


    GridForwardIterator<FunctionElement<dim,range,1>> elem(grid, 0);
    GridForwardIterator<FunctionElement<dim,range,1>> end(grid, IteratorState::pass_the_end);

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

