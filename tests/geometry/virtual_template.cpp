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
template class FunctionElement<2, 0, 1, 1> ;
template class GridForwardIterator<FunctionElement<2, 0, 1, 1>> ;


//class function {
//public:
//    template<int k>
//    virtual void reset() const;
//};

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



//using v = seq<Quadrature, 1, 2>::type;//boost::mpl::vector<Quadrature<1>,Quadrature<2>>;

template<int dim>
class Function : public GridElementHandler<dim>
{
    using v = typename seq<Quadrature, dim-1, dim>::type;
    using variant_1 = typename boost::make_variant_over<v>::type;
    using variant_type = boost::variant<
            Quadrature<1>,
            Quadrature<2> > ;

    using ElementAccessor = FunctionElement<dim, 0, 1, 1>;
    using ElementIterator = GridForwardIterator<ElementAccessor>;
    using parent_t = GridElementHandler<dim>;
public:
    Function(std::shared_ptr<const CartesianGrid<dim>> grid)
:
    GridElementHandler<dim>(grid)
    {}

    void init_cache(ElementIterator &elem, const variant_1& quad)
    {
        init_cache(elem.get_accessor(), quad);
    }

    void init_cache(ElementAccessor &elem, const variant_1& quad)
    {
        init_cache_impl.grid = this;
        init_cache_impl.elem = &elem;
        boost::apply_visitor(init_cache_impl, quad);
    }

    void fill_cache(ElementIterator &elem, const int j, const variant_type& quad)
    {
        fill_cache(elem.get_accessor(), j, quad);
    }

    void fill_cache(ElementAccessor &elem, const int j, const variant_type& quad)
    {
        fill_cache_impl.j = j;
        fill_cache_impl.grid = this;
        fill_cache_impl.elem = &elem;
        boost::apply_visitor(fill_cache_impl, quad);
    }

    virtual void reset(const NewValueFlags &flag, const variant_type& quad)
    {
        reset_impl.flag = flag;
        reset_impl.grid = this;
        boost::apply_visitor(reset_impl, quad);
    }

private:
    struct ResetDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T& quad)
        {
            grid->template reset<T::dim>(flag, quad);
        }

        NewValueFlags flag;
        parent_t *grid;
    };

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T& quad)
        {
            grid->template fill_cache<T::dim>(*elem, j);
        }

        int j;
        parent_t *grid;
        ElementAccessor *elem;
    };

    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T& quad)
        {
            grid->template init_cache<T::dim>(*elem);
        }

        parent_t *grid;
        ElementAccessor *elem;
    };

    ResetDispatcher reset_impl;

    FillCacheDispatcher fill_cache_impl;

    InitCacheDispatcher init_cache_impl;
};


int main() {
    const int dim = 2;
    auto grid = CartesianGrid<dim>::create(3);
    Function<2> x(grid);

    GridForwardIterator<FunctionElement<dim,0,1,1>> elem(grid, 0);
    //    GridForwardIterator<FunctionElement<dim,0,1,1>> end(grid, IteratorState::pass_the_end);
    //
    //    x.reset(NewValueFlags::none, QGauss<1>(2));
    x.reset(NewValueFlags::value, QGauss<2>(2));
    x.init_cache(elem, QGauss<2>(2));
    x.fill_cache(elem, 0, QGauss<2>(2));

}

