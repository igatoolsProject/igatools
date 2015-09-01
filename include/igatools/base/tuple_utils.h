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

#ifndef __TUPLE_UTILS_H_
#define __TUPLE_UTILS_H_

#include <igatools/base/config.h>
#include <igatools/geometry/unit_element.h>
#include <igatools/base/quadrature.h>

#include <tuple>

#include <boost/mpl/int.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/make_map.hpp>
#include <boost/fusion/include/for_each.hpp>

#include <memory>
IGA_NAMESPACE_OPEN

/**
 *
 * Returns a boost::fusion::map in which the keys are the <em>types</em>
 * <tt>Topology<I_min>,Topology<I_min+1>,Topology<I_min+2>,...</tt>
 * and the associated values are the <em>objects</em> of the type
 * <tt>DataIndexed<I_min>,DataIndexed<I_min+1>,DataIndexed<I_min+2>,...</tt>.
 *
 * The last index used is given by the sum of <tt>I_min</tt> with the size of
 * the index sequence used as input parameter of the function.
 *
 * It is used to infer (at compile time) the type of a boost::fusion::map when
 * the keys and values are two classes indexed with an integer value.
 *
 * @warning This function is implemented using some ''<em>black magic</em>''
 * template metaprogramming techniques.
 *
 * \tparam DataIndexed Template class representing the object that has to be indexed and used as
 * <tt>value type</tt> in the boost::fusion::map container. The template argument of the class must be
 * the integer used also for the indexing.
 * \tparam I_min Minimum index value
 */
template<template <int> class DataIndexed,int I_min,std::size_t... I>
auto
make_fusion_map_indexed_data(std::index_sequence<I...>)
{
  return boost::fusion::map<
         boost::fusion::pair<Topology<I+I_min>,DataIndexed<I+I_min> > ...>(
           boost::fusion::pair<Topology<I+I_min>,DataIndexed<I+I_min> >() ...);
}

/**
 * Alias for a boost::fusion::map container, in which the keys are the <em>types</em>
 * <tt>Topology<I_min>,Topology<I_min+1>,Topology<I_min+2>,...,Topology<I_min+N-1></tt>
 * and the associated values are the <em>objects</em> of the type
 * <tt>DataSameId<I_min>,DataSameId<I_min+1>,DataSameId<I_min+2>,...,DataSameId<I_min+N-1></tt>.
 *
 * * \tparam DataSameId Template class representing the object that has to be indexed and used as
 * <tt>value type</tt> in the boost::fusion::map container. The template argument of the class must be
 * the integer used also for the indexing.
 * \tparam I_min Minimum index value
 * \tparam N Number of objects in the boost::fusion::map container.
 *
 * @note The indices used run from <tt>I_min</tt> to <tt>I_min+N-1</tt>.
 *
 * @sa make_fusion_map_indexed_data
 */
template <template <int> class DataSameId,int Id_min,int N>
using DataVaryingId = decltype(make_fusion_map_indexed_data<DataSameId,Id_min>(std::make_index_sequence<N>()));



//TODO (pauletti, Aug 18, 2015): put this class in another file?
/**
 * List of Quadrature for the sub-elements having their topological dimension
 * ranging from <tt>dim-num_sub_elem</tt> to <tt>dim</tt>
 *
 * @note <tt>num_sub_elem</tt> is defined at configuration time in the main
 * CMakeLists.txt file.
 * @ingroup serializable
 */
template <int dim>
using QuadPtr = std::shared_ptr<const Quadrature<dim>>;

template<int dim>
class QuadList
  : public DataVaryingId<QuadPtr, (num_sub_elem <= dim ? dim - num_sub_elem : dim), (num_sub_elem <= dim ? num_sub_elem+1 : 1)>
{
public:
  template<int sdim>
  const QuadPtr<sdim> &get_quad() const
  {
    const auto &quad = boost::fusion::at_key<Topology<sdim>>(*this);
    Assert(quad != nullptr,ExcNullPtr());
    return quad;
  }

  template<int sdim>
  QuadPtr<sdim> &get_quad()
  {
    return boost::fusion::at_key<Topology<sdim>>(*this);
  }

private:
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class boost::serialization::access;

  template<class Archive>
  void
  serialize(Archive &ar, const unsigned int version)
  {
    using namespace boost::serialization;
    using namespace boost::fusion;
    for_each(*this,
             [&](auto & quad_same_topology_dim)
    {


      using PairType = typename std::remove_reference<decltype(quad_same_topology_dim)>::type;
      using SubDimType = typename PairType::first_type;
      std::string tag_name = "quad_" + std::to_string(SubDimType::value);
      auto non_const = std::const_pointer_cast<Quadrature<SubDimType::value>>(quad_same_topology_dim.second);

      ar &make_nvp(tag_name.c_str(), non_const);
      quad_same_topology_dim.second = non_const;
      Assert(quad_same_topology_dim.second != nullptr, ExcNullPtr());

    }
            );

  };
  ///@}
};

IGA_NAMESPACE_CLOSE



#endif
