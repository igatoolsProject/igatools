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


#ifndef VALUES_CACHE_H_
#define VALUES_CACHE_H_

#include <igatools/base/config.h>
#include <igatools/base/value_types.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/flags_handler.h>

#include <igatools/utils/value_vector.h>
#include <igatools/utils/value_table.h>

#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/any.hpp>

IGA_NAMESPACE_OPEN

/**
 * @ingroup cache
 */
template <class ValuesCache,int dim,int sdim>
using ValuesCacheAllSubElems = SafeSTLArray<ValuesCache,UnitElement<dim>::template num_elem<sdim>()>;

/**
 * @ingroup cache
 */
template<class ValuesCache, int dim, std::size_t... I>
auto
make_fusion_map_cache_all_sub_elems(std::index_sequence<I...>)
{
  return boost::fusion::map<boost::fusion::pair<
         Topology<(dim>I) ? dim-I : 0>,
         ValuesCacheAllSubElems<ValuesCache,dim,(dim>I) ? dim-I : 0> > ...>(
           boost::fusion::pair<Topology<(dim>I) ? dim-I : 0>,
           ValuesCacheAllSubElems<ValuesCache,dim,(dim>I) ? dim-I : 0> >() ...);
}


/**
 * @brief List of ValuesCache for the sub-elements having their topological dimension
 * ranging from <tt>dim</tt> to <tt>dim-num_sub_elem+1</tt>
 *
 * @note <tt>num_sub_elem</tt> is defined at configuration time in the main
 * CMakeLists.txt file.
 *
 * @ingroup cache
 */
template <class ValuesCache,int dim>
using DataCacheAllSubElems = decltype(make_fusion_map_cache_all_sub_elems<ValuesCache,dim>(
                                        std::make_index_sequence<dim+1>()));
//                                        std::make_index_sequence<(num_sub_elem <= dim ? num_sub_elem+1 : 1)>()));








/**
 * @ingroup cache
 */
namespace cacheutils
{

template<class FusionContainer>
void
print_caches(const FusionContainer &data, LogStream &out)
{
  boost::fusion::for_each(data,
                          [&](const auto & data_same_topology_dim)
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


#if 0
/**
 * @ingroup cache
 */
template<class CacheType>
ValueFlags
get_valid_flags_from_cache_type(const CacheType &cache)
{
  ValueFlags valid_flags = ValueFlags::none;

  boost::fusion::for_each(cache,
                          [&](const auto & type_and_status) -> void
  {
    using ValueType_Status = typename std::remove_reference<decltype(type_and_status)>::type;
    using ValueType = typename ValueType_Status::first_type;

    valid_flags |= ValueType::flag;
  } // end lambda function
                         );
  return valid_flags;

}
#endif

template <int sub_elem_dim, class FusionContainer>
const auto &
extract_sub_elements_data(const FusionContainer &data)
{
  return boost::fusion::at_key<Topology<sub_elem_dim>>(data);
}

template <int sub_elem_dim, class FusionContainer>
auto &
extract_sub_elements_data(FusionContainer &data)
{
  return boost::fusion::at_key<Topology<sub_elem_dim>>(data);
}




} // end namespace cacheutils


/**
 * @brief Numerical values (ValueTable or a ValueVector) with an associated FlagStatus.
 * It is used to store the values of a given quantity
 * in the element's cache (for a given dimension and for a given sub-element).
 *
 * @tparam DataType This should be a ValueTable or a ValueVector.
 *
 * @ingroup cache
 */
template < class DataType >
class DataWithFlagStatus : public DataType
{
  using self_t = DataWithFlagStatus<DataType>;

public:
  using DataType::DataType;


  /**
   * Fills the internal values with @p data and
   * sets the associated <b>filled</b> status to <b>true</b>.
   */
  void fill(const DataType &data)
  {
    Assert(status_.fill_,ExcMessage("The data is not intended to be filled."));
    if (this != &data)
    {
      DataType::operator=(data);
      status_.filled_ = true;
    }
  }

  /**
   * Fills the internal values with @p data and
   * sets the associated <b>filled</b> status to <b>true</b>.
   */
  void fill(DataType &&data)
  {
    Assert(status_.fill_,ExcMessage("The data is not intended to be filled."));
    DataType::operator=(std::move(data));
    status_.filled_ = true;
  }

  bool status_fill() const
  {
    return status_.fill_;
  };

  void set_status_fill(const bool fill_status)
  {
    status_.fill_ = fill_status;
  };

  bool status_filled() const
  {
    return status_.filled_;
  };


  void set_status_filled(const bool filled_status)
  {
    status_.filled_ = filled_status;
  };
//*/

  const FlagStatus &get_status() const
  {
    return status_;
  }

  /**
   * Returns an estimate of the memory used to define the object.
   */
  auto memory_consumption() const
  {
    return DataType::memory_consumption() + sizeof(status_) ;
  }

private:
  FlagStatus status_;
};

/**
 *
 * @ingroup cache
 */
template<int dim,class CacheType>
class ValuesCache : public CacheStatus
{
public:
  static constexpr int get_dim()
  {
    return dim;
  }

protected:
  CacheType values_;

public:
  CacheType &get_values()
  {
    return values_;
  }

  /**
   * Function for printing information of the cache.
   * Its use is intended for testing and debugging purposes.
   */
  void print_info(LogStream &out) const
  {
    out.begin_item("Memory consumption: " + std::to_string(this->memory_consumption()) + " bytes");
    out.end_item();

    boost::fusion::for_each(values_,
                            [&](const auto & type_and_value) -> void
    {
      using ValueType_ValueContainer = typename std::remove_reference<decltype(type_and_value)>::type;
      using ValueType = typename ValueType_ValueContainer::first_type;
      const auto &value = type_and_value.second;

      out.begin_item(ValueType::name + ": (memory consumption: " + std::to_string(value.memory_consumption()) + " bytes)");

      out.begin_item("Fill flags:");
      value.get_status().print_info(out);
      out.end_item();


      if (value.status_filled())
      {
        out.begin_item("Data:");
        value.print_info(out);
        out.end_item();
      }
      out.end_item();

    } // end lambda function
                           );
  }

  /**
   * Returns an estimate of the memory used to define the object.
   */
  auto memory_consumption() const
  {
    std::size_t memory_consumption = 0;
    boost::fusion::for_each(values_,
                            [&](const auto & type_and_value) -> void
    {
      memory_consumption += type_and_value.second.memory_consumption();
    } // end lambda function
                           );
    return memory_consumption;
  }


  /**
   * Returns the data associated with the given <tt>ValueType</tt> (non-const version).
   */
  template<class ValueType>
  auto &get_data()
  {
    return boost::fusion::at_key<ValueType>(values_);
  }

  /**
   * Returns the data associated with the given <tt>ValueType</tt> (const version).
   *
   * @note In Debug mode, a check is performed in order to check that the FlagStatus of the
   * requested data is marked as "filled", otherwise, an assertion will be raised.
   */
  template<class ValueType>
  const auto &get_data() const
  {
    const auto &data = boost::fusion::at_key<ValueType>(values_);
    Assert(data.status_filled(),
           ExcMessage("The cache for \"" + ValueType::name + "\" is not filled."));

    return data;
  }


  /**
   * @name Functions used to query or modify the Flag status for a given ValueType
   */
  ///@{
  /** Returns true if the quantity associated to @p ValueType must be filled. */
  template<class ValueType>
  bool status_fill() const
  {
    return boost::fusion::at_key<ValueType>(values_).status_fill();
  }

  /** Returns true if the quantity associated to @p ValueType is filled. */
  template<class ValueType>
  bool status_filled() const
  {
    return boost::fusion::at_key<ValueType>(values_).status_filled();
  }


  /** Returns true if the nothing must be filled. */
  bool fill_none() const
  {
    const bool fill_someone = boost::fusion::any(values_,
                                                 [](const auto & type_and_data) -> bool
    {
      return type_and_data.second.status_fill() == true;
    } // end lambda function
                                                );
    return !fill_someone;
  }
  ///@}

};



/**
 * @brief This class is used to wrap the cache for values associated basis fucnctions
 * and to evaluation points.
 *
 *
 * @ingroup cache
 */
template<int dim, class CacheType>
class BasisValuesCache : public ValuesCache<dim,CacheType>
{

public:
  /**
   * Allocate space for the values and derivatives
   * of the element basis functions at quadrature points
   * as specify by the @p flags.
   */
  void resize(const typename space_element::Flags &flags,
              const Size n_points,
              const Size n_basis)
  {
    Assert(n_points >= 0, ExcLowerRange(n_points,1));
    Assert(n_basis > 0, ExcLowerRange(n_basis,1));

    boost::fusion::for_each(this->values_,
                            [&](auto & type_and_value) -> void
    {
      using ValueType_ValueContainer = typename std::remove_reference<decltype(type_and_value)>::type;
      using ValueType = typename ValueType_ValueContainer::first_type;
      auto &value = type_and_value.second;

      if (contains(flags,ValueType::flag))
      {
        value.set_status_fill(true);

        if (value.get_num_points() != n_points ||
        value.get_num_functions() != n_basis)
        {
          value.resize(n_basis,n_points);
        }
        value.zero();
      }
      else
      {
        value.clear();
        value.set_status_fill(false);
      }
      value.set_status_filled(false);
    } // end lambda function
                           );

    this->set_initialized(true);
  }
};

/**
 * @brief This class is used to wrap the cache for values associated to evaluation points.
 *
 *
 * @ingroup cache
 */
template<int dim, class CacheType>
class PointValuesCache : public ValuesCache<dim,CacheType>
{

public:

  /**
   * Allocate space for the cache values
   * of the element at evaluation points
   * as specify by the @p flags.
   */
  template<class Flags>
  void resize(const Flags &flags, const Size n_points)
  {
    Assert(n_points >= 0, ExcLowerRange(n_points,1));

    boost::fusion::for_each(this->values_,
                            [&](auto & type_and_value) -> void
    {
      using ValueType_ValueContainer =
      typename std::remove_reference<decltype(type_and_value)>::type;
      using ValueType = typename ValueType_ValueContainer::first_type;
      auto &value = type_and_value.second;

      if (contains(flags, ValueType::flag))
      {
        value.set_status_fill(true);

        if (value.get_num_points() != n_points)
          value.resize(n_points);

        value.zero();
      }
      else
      {
        value.clear();
        value.set_status_fill(false);
      }
      value.set_status_filled(false);

    } // end lambda function
                           );

    this->set_initialized(true);
  }
};




/**
 *
 *
 * @ingroup cache
 *
 */
template <class SubElemCache>
class AllSubElementsCache
{
public:
  AllSubElementsCache() = default;

  AllSubElementsCache(const AllSubElementsCache &in) = default;

  AllSubElementsCache(AllSubElementsCache &&in) = default;

  ~AllSubElementsCache() = default;


  AllSubElementsCache &operator=(const AllSubElementsCache &in) = delete;

  AllSubElementsCache &operator=(AllSubElementsCache &&in) = delete;



  void print_info(LogStream &out) const
  {
    out.begin_item("Cache for all sub-elements in all dimensions: (memory consumption: " + std::to_string(this->memory_consumption()) + " bytes)");

    boost::fusion::for_each(cache_all_sub_elems_,
                            [&](const auto & data_same_topology_dim)
    {
      using PairType = typename std::remove_reference<decltype(data_same_topology_dim)>::type;
      using SubDimType = typename PairType::first_type;

      const auto &data_same_subdim = data_same_topology_dim.second;

      out.begin_item("Cache for sub-element(s) with dimension: " + std::to_string(SubDimType::value) +
                     "  (memory consumption : " +
                     std::to_string(data_same_subdim.size() * data_same_subdim[0].memory_consumption()) + " bytes.");
      int s_id =0;
      for (const auto &data_same_topology_id : data_same_subdim)
      {
        out.begin_item("Sub-element id: " + std::to_string(s_id++));
        data_same_topology_id.print_info(out);
        out.end_item();
      }
      out.end_item();
    }
                           );


    out.end_item();
  }

  template <int sub_elem_dim>
  SubElemCache &
  get_sub_elem_cache(const int sub_elem_id)
  {
    return cacheutils::extract_sub_elements_data<sub_elem_dim>(cache_all_sub_elems_)[sub_elem_id];
  }

  template <int sub_elem_dim>
  const SubElemCache &
  get_sub_elem_cache(const int sub_elem_id) const
  {
    const auto &cache = cacheutils::extract_sub_elements_data<sub_elem_dim>(cache_all_sub_elems_)[sub_elem_id];
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    return cache;
  }

#if 0
  ValueFlags get_valid_flags() const
  {
    return this->template get_sub_elem_cache<SubElemCache::get_dim()>(0).get_valid_flags();
  }
#endif

  /**
   * Returns an estimate of the memory used to define the object.
   */
  auto memory_consumption() const
  {
    std::size_t memory_consumption = 0;
    boost::fusion::for_each(cache_all_sub_elems_,
                            [&](const auto & type_and_value) -> void
    {
      for (const auto &cache : type_and_value.second)
        memory_consumption += cache.memory_consumption();
    } // end lambda function
                           );
    return memory_consumption;
  }


  /**
   * Cache for all sub-elements.
   */
  DataCacheAllSubElems<SubElemCache, SubElemCache::get_dim()> cache_all_sub_elems_;
};


IGA_NAMESPACE_CLOSE



#endif // #ifndef VALUES_CACHE_H_

