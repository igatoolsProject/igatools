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

#include <boost/mpl/int.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/make_map.hpp>
#include <boost/fusion/include/for_each.hpp>


IGA_NAMESPACE_OPEN

template<template<int> class Q, std::size_t... I>
auto
tuple_of_quads(std::index_sequence<I...>)
{
    return boost::fusion::make_vector(Q<I>() ...);
}

template<int dim, template<int> class Q>
using TupleList = decltype(tuple_of_quads<Q>(std::make_index_sequence<dim+1>()));

template<int dim>
using EvalPtsList = TupleList<dim, Quadrature>;


template<class ValuesCache, int dim, std::size_t... I>
auto
tuple_of_caches(
    std::index_sequence<I...>,
    const Quadrature<dim> &q,
    const ValuesCache &)
{
    using Topology = UnitElement<dim>;
    return boost::fusion::map<boost::fusion::pair<boost::mpl::int_<(dim>I) ? dim-I : 0>,std::array<ValuesCache,Topology::template num_elem<(dim>I) ? dim-I : 0>()>> ...>(
               boost::fusion::pair<boost::mpl::int_<(dim>I) ? dim-I : 0>,std::array<ValuesCache,Topology::template num_elem<(dim>I) ? dim-I : 0>()>>() ...);
}


template<class ValuesCache, int dim>
using CacheList = decltype(tuple_of_caches(
                               std::make_index_sequence<(num_sub_elem <= dim ? num_sub_elem+1 : 1)>(),
                               Quadrature<dim>(),
                               ValuesCache()));





namespace cacheutils
{

template<class FusionContainer>
void print_caches(const FusionContainer &data, LogStream &out)
{
    boost::fusion::for_each(data,
                            [&out](const auto & data_same_topology_dim)
    {
        using PairType = typename std::remove_reference<decltype(data_same_topology_dim)>::type;
        using SubDimType = typename PairType::first_type;
        out.begin_item("Cache for sub-element(s) with dimension: " + std::to_string(SubDimType::value));
        int s_id =0;
        for (const auto &data_same_topology_id : data_same_topology_dim.second)
        {
            out.begin_item("Sub-element id: " + std::to_string(s_id++));
            data_same_topology_id.print_info(out);
            out.end_item();
        }
        out.end_item();
    }
                           );
}

template <int sub_elem_dim, class FusionContainer>
const auto &
extract_sub_elements_data(const FusionContainer &data)
{
    return boost::fusion::at_key<boost::mpl::int_<sub_elem_dim>>(data);
}

template <int sub_elem_dim, class FusionContainer>
auto &
extract_sub_elements_data(FusionContainer &data)
{
    return boost::fusion::at_key<boost::mpl::int_<sub_elem_dim>>(data);
}

}; // end namespace cacheutils


IGA_NAMESPACE_CLOSE

#endif
