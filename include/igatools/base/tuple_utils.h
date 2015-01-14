//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

#ifndef __IGA_TUPLE_UTILS_H_
#define __IGA_TUPLE_UTILS_H_

#include <igatools/base/config.h>
#include <igatools/geometry/unit_element.h>
#include <igatools/base/quadrature.h>
#include <tuple>

IGA_NAMESPACE_OPEN

template<template<int> class Q, std::size_t... I>
auto tuple_of_quads(std::index_sequence<I...>)
-> decltype(std::make_tuple(Q<I>() ...))
{
    return std::make_tuple(Q<I>() ...);
}

template<int dim, template<int> class Q>
using TupleList = decltype(tuple_of_quads<Q>(std::make_index_sequence<dim+1>()));

template<class ValuesCache, int dim, std::size_t... I>
auto tuple_of_caches(std::index_sequence<I...>, const QuadratureTensorProduct<dim> &q, const ValuesCache &)
-> decltype(std::make_tuple(std::array<ValuesCache,
                            UnitElement<dim>::template num_elem<I>()>() ...))
{
    return std::make_tuple(std::array<ValuesCache,
                           UnitElement<dim>::template num_elem<I>()>() ...);
}


template<class ValuesCache, int dim>
using CacheList = decltype(tuple_of_caches(std::make_index_sequence<dim+1>(),
                                           QuadratureTensorProduct<dim>(),
                                           ValuesCache()));

template<class Func, class Tuple, std::size_t N, std::size_t Min>
struct TupleFunc
{
    static void apply_func(Func &F, const Tuple &t)
    {
        TupleFunc<Func,Tuple, N-1, Min>::apply_func(F,t);
        if (N>Min)
            F.func(std::get<N-1>(t));
    }
};

template<class Func, class Tuple, std::size_t N>
struct TupleFunc<Func, Tuple, N, N>
{
    static void apply_func(Func &F, const Tuple &t)
    {
        F.func(std::get<N>(t));
    }
};


template<class Func, class Args1, class Args2, class Tuple, std::size_t N, std::size_t Min>
struct TupleFunc1
{
    static void apply_func(Func &F, const Args1 &flag, const Args2 &quad, Tuple &t)
    {
        TupleFunc1<Func, Args1, Args2, Tuple, N-1, Min>::apply_func(F, flag, quad, t);
        if (N>Min)
        {
            auto &val_cache = std::get<N-1>(t);
            int j=0;
            for (auto &s_cache : val_cache)
            {
                F.func(s_cache, flag, quad.template collapse_to_sub_element<N-1>(j));
                ++j;
            }
        }
    }
};

template<class Func, class Args1, class Args2, class Tuple, std::size_t N>
struct TupleFunc1<Func, Args1, Args2, Tuple, N, N>
{
    static void apply_func(Func &F, const Args1 &flag, const Args2 &quad, Tuple &t)
    {
        auto &val_cache = std::get<N>(t);
        int j=0;
        for (auto &s_cache : val_cache)
        {
            F.func(s_cache, flag, quad.template collapse_to_sub_element<N>(j));
            ++j;
        }
    }
};



namespace cacheutils
{
struct PrintCacheFunc
{
    PrintCacheFunc(LogStream &out1)
        :out(out1)
    {}

    void func(const auto &c)
    {
        for (auto &e : c)
            e.print_info(out);
        out << std::endl;
    }
    LogStream &out;
};

template<class... Args>
void print_caches(const std::tuple<Args...> &t, LogStream &out)
{
    PrintCacheFunc f(out);
    TupleFunc<PrintCacheFunc, decltype(t), sizeof...(Args), 0>::apply_func(f,t);
}
};


IGA_NAMESPACE_CLOSE

#endif
