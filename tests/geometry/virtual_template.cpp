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
 *  Test of
 *  author: pauletti
 *  date: 2014-08-18
 *
 */

#include "../tests.h"
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid_element_handler.h>
//#include <igatools/base/function_element.h>
#include <boost/variant.hpp>

#include <igatools/base/new_function-template.h>
#include <igatools/base/function_element.cpp>
#include <igatools/geometry/grid_forward_iterator.cpp>

#include <boost/mpl/vector.hpp>


template class NewFunction<2, 0, 1, 1> ;
//template class FunctionElement<2, 0, 1, 1> ;
//template class GridForwardIterator<FunctionElement<2, 0, 1, 1>> ;


//class function {
//public:
//    template<int k>
//    virtual void reset() const;
//};


template<int k_>
struct Int
{
    static const int k = k_;
};

template<template<int> class Q, int start, std::size_t N>
struct seq;

template<template<int> class Q, int start>
struct seq<Q, start, start>
{
    using type = boost::mpl::vector<Q<start>>;
};



template<template<int> class Q, int start, std::size_t N>
struct seq
{
    using v1 = typename seq<Q, start, N-1>::type;
    using type = typename boost::mpl::push_back<v1, Q<N>>::type;
};



template<int dim>
class Function : public GridElementHandler<dim>
{
    using v1 = typename seq<Quadrature, dim-1, dim>::type;
    using variant_1 = typename boost::make_variant_over<v1>::type;

    using v2 = typename seq<Int, dim-1, dim>::type;
    using variant_2 = typename boost::make_variant_over<v2>::type;

    using ElementAccessor = FunctionElement<dim, 0, 1, 1>;
    using ElementIterator = GridForwardIterator<ElementAccessor>;
    using parent_t = GridElementHandler<dim>;
public:
    Function(std::shared_ptr<const CartesianGrid<dim>> grid)
        :
        GridElementHandler<dim>(grid)
    {}

    void init_cache(ElementIterator &elem, const variant_2 &k)
    {
        init_cache(elem.get_accessor(), k);
    }

    void init_cache(ElementAccessor &elem, const variant_2 &k)
    {
        init_cache_impl.grid_handler = this;
        init_cache_impl.elem = &elem;
        boost::apply_visitor(init_cache_impl, k);
    }

    void fill_cache(ElementIterator &elem, const int j, const variant_2 &k)
    {
        fill_cache(elem.get_accessor(), j, k);
    }

    void fill_cache(ElementAccessor &elem, const int j, const variant_2 &k)
    {
        fill_cache_impl.j = j;
        fill_cache_impl.grid_handler = this;
        fill_cache_impl.elem = &elem;
        boost::apply_visitor(fill_cache_impl, k);
    }

    virtual void reset(const NewValueFlags &flag, const variant_1 &quad)
    {
        reset_impl.flag = flag;
        reset_impl.grid_handler = this;
        boost::apply_visitor(reset_impl, quad);
    }

private:
    struct ResetDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            grid_handler->template reset<T::dim>(flag, quad);
        }

        NewValueFlags flag;
        parent_t *grid_handler;
    };

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            grid_handler->template fill_cache<T::k>(*elem, j);
        }

        int j;
        parent_t *grid_handler;
        ElementAccessor *elem;
    };

    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            grid_handler->template init_cache<T::k>(*elem);
        }

        parent_t *grid_handler;
        ElementAccessor *elem;
    };

    ResetDispatcher reset_impl;
    FillCacheDispatcher fill_cache_impl;
    InitCacheDispatcher init_cache_impl;
};



template <int dim, int k>
void test()
{
    auto grid = CartesianGrid<dim>::create(3);
    Function<dim> x(grid);

    GridForwardIterator<FunctionElement<dim,0,1,1>> elem(grid, 0);
    GridForwardIterator<FunctionElement<dim,0,1,1>> end(grid, IteratorState::pass_the_end);

    x.reset(NewValueFlags::value, QGauss<k>(2));
    x.init_cache(elem, Int<k>());
    x.fill_cache(elem, 0, Int<k>());

    //elem->get_points().print_info(out);
}

int main()
{
    test<2,2>();
    test<2,1>();
}

