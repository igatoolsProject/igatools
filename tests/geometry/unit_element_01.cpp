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

/*
 *  Test of unit element
 *  author: pauletti
 *  date: 2014-08-18
 *
 */

#include "../tests.h"

#include <igatools/geometry/unit_element.h>
#include <igatools/base/quadrature.h>

template<template<int> class Q, std::size_t... I>
auto tuple_of_quads(std::index_sequence<I...>)
-> decltype(std::make_tuple(Q<I>() ...))
{
    return std::make_tuple(Q<I>() ...);
}

template<int dim, template<int> class Q>
using TupleList = decltype(tuple_of_quads<Q>(std::make_index_sequence<dim+1>()));

template<int dim>
using QuadList = TupleList<dim, QuadratureTensorProduct>;



template<class ValuesCache, int dim, std::size_t... I>
auto tuple_of_caches(std::index_sequence<I...>, const QuadratureTensorProduct<dim> &q, const ValuesCache &)
-> decltype(std::make_tuple(std::array<ValuesCache,
                            UnitElement<dim>::template num_elem<dim-I>()>() ...))
{
    return std::make_tuple(std::array<ValuesCache,
                           UnitElement<dim>::template num_elem<dim-I>()>() ...);
}


template<int dim, int n_sub_elem>
using CacheList = decltype(tuple_of_caches(std::make_index_sequence<n_sub_elem+1>(),
                                           QuadratureTensorProduct<dim>(),
                                           3));


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




struct PrintQuadFunc
{
    PrintQuadFunc(LogStream &out1)
        :out(out1)
    {}
    template<std::size_t... I>
    void func(const auto &q)
    {
        q.print_info(out);
        out << endl;
    }
    LogStream &out;
};

template<class... Args>
void print_quads(const std::tuple<Args...> &t, LogStream &out1)
{
    PrintQuadFunc f(out1);
    TupleFunc<PrintQuadFunc, decltype(t), sizeof...(Args), 2>::apply_func(f,t);
}




int main()
{



    const int dim=3;
    QuadList<dim>  list_of_quad;
    print_quads(list_of_quad, out);

#if 1

    CacheList<3,1> list_values;
    auto &cache0 = std::get<0>(list_values);
    auto &cache1 = std::get<1>(list_values);

    for (auto &c : cache0)
        out << c << endl;

    for (auto &c : cache1)
        out << c << endl;

    UnitElement<2>::num_elem<1>();
    UnitElement<2>::get_elem<1>(0);
#endif



    return 0;
}
