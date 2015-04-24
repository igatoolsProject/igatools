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

#include <igatools/base/function.h>

#include <igatools/base/quadrature.h>

#include <igatools/utils/value_vector.h>
#include <igatools/utils/value_table.h>



#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/at_key.hpp>

IGA_NAMESPACE_OPEN

#if 0
template<int dim, int sub_elem_dim, class CacheSingleElem>
class CacheSubElems : std::array<CacheSingleElem,UnitElement<dim>::template num_elem<sub_elem_dim>()>
{
    using self_t = std::array<CacheSingleElem,UnitElement<dim>::template num_elem<sub_elem_dim>()>;
public:

    using self_t::self_t;
};
#endif

template<class ValuesCache, int dim, std::size_t... I>
auto
tuple_of_caches(
    std::index_sequence<I...>,
    const Quadrature<dim> &q,
    const ValuesCache &)
{
    return boost::fusion::map<boost::fusion::pair<Topology<(dim>I) ? dim-I : 0>,std::array<ValuesCache,UnitElement<dim>::template num_elem<(dim>I) ? dim-I : 0>()>> ...>(
               boost::fusion::pair<Topology<(dim>I) ? dim-I : 0>,std::array<ValuesCache,UnitElement<dim>::template num_elem<(dim>I) ? dim-I : 0>()>>() ...);
}


template<class ValuesCache, int dim>
using CacheList = decltype(tuple_of_caches(
                               std::make_index_sequence<(num_sub_elem <= dim ? num_sub_elem+1 : 1)>(),
                               Quadrature<dim>(),
                               ValuesCache()));







template < class DataType >
class DataWithFlagStatus : public DataType
{
    using self_t = DataWithFlagStatus<DataType>;

public:
    using DataType::DataType;


    DataWithFlagStatus &operator=(const DataType &data)
    {
        if (this != &data)
        {
            DataType::operator=(data);
            status_.filled_ = true;
        }
        return (*this);
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

    const FlagStatus &get_status() const
    {
        return status_;
    }
//*/
    /*
        void set_status(const FlagStatus &status)
        {
            status_ = status;
        }
    //*/

private:
//  DataType data_;
    FlagStatus status_;
};


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
    void print_info(LogStream &out) const
    {
        boost::fusion::for_each(values_,
                                [&](const auto & type_and_value) -> void
        {
            using ValueType_ValueContainer = typename std::remove_reference<decltype(type_and_value)>::type;
            using ValueType = typename ValueType_ValueContainer::first_type;
            out.begin_item(ValueType::name + ":");

            out.begin_item("Fill flags:");
            type_and_value.second.get_status().print_info(out);
            out.end_item();

            if (type_and_value.second.status_filled())
            {
                out.begin_item("Data:");
                type_and_value.second.print_info(out);
                out.end_item();
            }
            out.end_item();
        } // end lambda function
                               );
    }


    template<class ValueType>
    auto &get_data()
    {
        return boost::fusion::at_key<ValueType>(values_);
    }

    template<class ValueType>
    const auto &get_data() const
    {
        const auto &data = boost::fusion::at_key<ValueType>(values_);
        //TODO (martinelli, Apr 03,2015): uncomment this assertion
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

    /** Sets the filled @p status the quantity associated to @p ValueType. */
    template<class ValueType>
    void set_status_filled(const bool status)
    {
        this->template get_data<ValueType>().set_status_filled(status);
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



    /**
     * Returns the flags that are valid to be used with this class.
     *
     * @note The valid flags are defined to be the ones that can be inferred from the ValueType(s)
     * used as key of the boost::fusion::map in CacheType.
     */
    ValueFlags get_valid_flags() const
    {
        ValueFlags valid_flags = ValueFlags::none;

        boost::fusion::for_each(values_,
                                [&](const auto & type_and_status) -> void
        {
            using ValueType_Status = typename std::remove_reference<decltype(type_and_status)>::type;
            using ValueType = typename ValueType_Status::first_type;

            valid_flags |= ValueType::flag;
        } // end lambda function
                               );
        return valid_flags;
    }

};




template<int dim, class CacheType>
class BasisValuesCache : public ValuesCache<dim,CacheType>
{

public:
    /**
     * Allocate space for the values and derivatives
     * of the element basis functions at quadrature points
     * as specify by the flag
     */
    void resize(const ValueFlags &flags,
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


template<int dim, class CacheType>
class FuncValuesCache : public ValuesCache<dim,CacheType>
{

public:

    /**
     * Allocate space for the fucntion values and derivatives
     * of the element at quadrature points
     * as specify by the flag
     */
    void resize(const ValueFlags &flags,
                const Size n_points)
    {
//        this->flags_handler_ = flags_handler;

        Assert(n_points >= 0, ExcLowerRange(n_points,1));

        boost::fusion::for_each(this->values_,
                                [&](auto & type_and_value) -> void
        {
            using ValueType_ValueContainer = typename std::remove_reference<decltype(type_and_value)>::type;
            using ValueType = typename ValueType_ValueContainer::first_type;
            auto &value = type_and_value.second;

            if (contains(flags,ValueType::flag))
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





template <class SubElemCache>
class LocalCache
{
public:
    LocalCache() = default;

    LocalCache(const LocalCache &in) = default;

    LocalCache(LocalCache &&in) = default;

    ~LocalCache() = default;


    LocalCache &operator=(const LocalCache &in) = delete;

    LocalCache &operator=(LocalCache &&in) = delete;

    void print_info(LogStream &out) const
    {
        cacheutils::print_caches(cache_all_sub_elems_, out);
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


    ValueFlags get_valid_flags() const
    {
        return this->template get_sub_elem_cache<SubElemCache::get_dim()>(0).get_valid_flags();
    }

    /**
     * Cache for all sub-elements.
     */
    CacheList<SubElemCache, SubElemCache::get_dim()> cache_all_sub_elems_;

};


IGA_NAMESPACE_CLOSE



#endif // #ifndef VALUES_CACHE_H_

