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

template<int dim,
         class CacheType,
         class FlagsType>
class ValuesCache : public CacheStatus
{
public:

    static constexpr int get_dim()
    {
        return dim;
    }



    FlagsType flags_handler_;

protected:


    CacheType values_;



public:
    void print_info(LogStream &out) const
    {
        out.begin_item("Fill flags:");
        flags_handler_.print_info(out);
        out.end_item();

        boost::fusion::for_each(values_,
                                [&](const auto & type_and_value) -> void
        {
            using ValueType_ValueContainer = typename std::remove_reference<decltype(type_and_value)>::type;
            using ValueType = typename ValueType_ValueContainer::first_type;
            if (flags_handler_.template filled<ValueType>())
            {
                out.begin_item(ValueType::name + "s:");
                type_and_value.second.print_info(out);
                out.end_item();
            }
        } // end lambda function
                               );
    }


    template<class ValueType>
    auto &get_der()
    {
        return boost::fusion::at_key<ValueType>(values_);
    }

    template<class ValueType>
    const auto &get_der() const
    {
        //TODO (martinelli, Apr 03,2015): uncomment this assertion
//        Assert(flags_handler_.template filled<ValueType>(),
//               ExcMessage("The cache for \"" + ValueType::name + "\" is not filled."));

        return boost::fusion::at_key<ValueType>(values_);
    }

};




template<int dim, class CacheType, class FlagsType>
class BasisValuesCache : public ValuesCache<dim,CacheType,FlagsType>
{

public:
    /**
     * Allocate space for the values and derivatives
     * of the element basis functions at quadrature points
     * as specify by the flag
     */
    void resize(const FlagsType &flags_handler,
                const Size n_points,
                const Size n_basis)
    {
        this->flags_handler_ = flags_handler;

        Assert(n_points >= 0, ExcLowerRange(n_points,1));
        Assert(n_basis > 0, ExcLowerRange(n_basis,1));

        boost::fusion::for_each(this->values_,
                                [&](auto & type_and_value) -> void
        {
            using ValueType_ValueContainer = typename std::remove_reference<decltype(type_and_value)>::type;
            using ValueType = typename ValueType_ValueContainer::first_type;
            auto &value = type_and_value.second;

            if (this->flags_handler_.template fill<ValueType>())
            {
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
                this->flags_handler_.template set_filled<ValueType>(false);
            }
        } // end lambda function
                               );

        this->set_initialized(true);
    }

};


template<int dim, class CacheType, class FlagsType>
class FuncValuesCache : public ValuesCache<dim,CacheType,FlagsType>
{

public:

    /**
     * Allocate space for the fucntion values and derivatives
     * of the element at quadrature points
     * as specify by the flag
     */
    void resize(const FlagsType &flags_handler,
                const Size n_points)
    {
        this->flags_handler_ = flags_handler;

        Assert(n_points >= 0, ExcLowerRange(n_points,1));

        boost::fusion::for_each(this->values_,
                                [&](auto & type_and_value) -> void
        {
            using ValueType_ValueContainer = typename std::remove_reference<decltype(type_and_value)>::type;
            using ValueType = typename ValueType_ValueContainer::first_type;
            auto &value = type_and_value.second;

            if (this->flags_handler_.template fill<ValueType>())
            {
                if (value.get_num_points() != n_points)
                    value.resize(n_points);

                value.zero();
            }
            else
            {
                value.clear();
                this->flags_handler_.template set_filled<ValueType>(false);
            }
        } // end lambda function
                               );

        this->set_initialized(true);
    }
};





template <class Cache>
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
        cacheutils::print_caches(values_, out);
    }

    template <int topology_dim>
    Cache &
    get_value_cache(const int topology_id)
    {
        return cacheutils::extract_sub_elements_data<topology_dim>(values_)[topology_id];
    }

    template <int topology_dim>
    const Cache &
    get_value_cache(const int topology_id) const
    {
        const auto &cache = cacheutils::extract_sub_elements_data<topology_dim>(values_)[topology_id];
        Assert(cache.is_filled() == true, ExcCacheNotFilled());
        return cache;
    }

    CacheList<Cache, Cache::get_dim()> values_;

};


IGA_NAMESPACE_CLOSE



#endif // #ifndef VALUES_CACHE_H_

