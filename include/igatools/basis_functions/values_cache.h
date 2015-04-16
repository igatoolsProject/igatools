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



#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/container/map/map_fwd.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/sequence/intrinsic/at_key.hpp>
#include <boost/fusion/include/at_key.hpp>

IGA_NAMESPACE_OPEN


template<int dim, int codim, int range, int rank, template<class ContainedType> class ContainerType>
class ValuesCache : public CacheStatus
{
    using Func = Function<dim,codim,range,rank>;

    using Value = typename Func::Value;

    template <int der_order>
    using Derivative = typename Func::template Derivative<der_order>;

    using Div = typename Func::Div;

public:

    static constexpr int get_dim()
    {
        return dim;
    }


    FunctionFlags flags_handler_;

    using CacheType = boost::fusion::map<
                      boost::fusion::pair<     _Value,ContainerType<Value>>,
                      boost::fusion::pair<  _Gradient,ContainerType<Derivative<1>>>,
                      boost::fusion::pair<   _Hessian,ContainerType<Derivative<2>>>,
                      boost::fusion::pair<_Divergence,ContainerType<Div>>
                      >;

protected:

    CacheType values_;



    struct ValuesCachePrinter
    {
        ValuesCachePrinter(const FunctionFlags &flags_handler,
                           LogStream &out)
            :
            flags_handler_ {flags_handler},
                       out_ {out}
        {}

        template <class pair_ValueType_ValueContainer>
        void operator()(const pair_ValueType_ValueContainer &x) const
        {
            using ValueType = typename pair_ValueType_ValueContainer::first_type;
            if (this->flags_handler_.template filled<ValueType>())
            {
                this->out_.begin_item(ValueType::name + "s:");
                x.second.print_info(this->out_);
                this->out_.end_item();
            }
        }


    private:
        const FunctionFlags &flags_handler_;

        LogStream &out_;
    };

public:
    void print_info(LogStream &out) const
    {
        out.begin_item("Fill flags:");
        flags_handler_.print_info(out);
        out.end_item();

        ValuesCachePrinter printer(flags_handler_,out);
        // apply the functor printer to each pair ValueType/Container in values_
        boost::fusion::for_each(values_,printer);
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
//        Assert(flags_handler_.filled<ValueType>(),
//               ExcMessage("The cache for " + ValueType::name + " is not filled."));

        return boost::fusion::at_key<ValueType>(values_);
    }

};




template<int dim, int codim, int range, int rank>
class BasisValuesCache : public ValuesCache<dim,codim,range,rank,ValueTable>
{
    using parent_t = ValuesCache<dim,codim,range,rank,ValueTable>;

private:
    struct BasisValuesCacheResizer
    {
        BasisValuesCacheResizer(
            FunctionFlags &flags_handler,
            const int n_basis,
            const int n_points)
            :
            flags_handler_ {flags_handler},
                       n_basis_ {n_basis},
        n_points_ {n_points}
        {
            Assert(n_points >= 0, ExcLowerRange(n_points,1));
            Assert(n_basis > 0, ExcLowerRange(n_basis,1));
        }

        template <class pair_ValueType_ValueContainer>
        void operator()(pair_ValueType_ValueContainer &x) const
        {
            using ValueType = typename pair_ValueType_ValueContainer::first_type;
            auto &value = x.second;

            if (this->flags_handler_.template fill<ValueType>())
            {
                if (value.get_num_points() != n_points_ ||
                    value.get_num_functions() != n_basis_)
                {
                    value.resize(n_basis_,n_points_);
                }
                value.zero();
            }
            else
            {
                value.clear();
                this->flags_handler_.template set_filled<ValueType>(false);
            }
        }


    private:
        FunctionFlags &flags_handler_;

        int n_basis_;
        int n_points_;
    };


public:
    /**
     * Allocate space for the values and derivatives
     * of the element basis functions at quadrature points
     * as specify by the flag
     */
    void resize(const FunctionFlags &flags_handler,
                const Size n_points,
                const Size n_basis)
    {
        this->flags_handler_ = flags_handler;

        BasisValuesCacheResizer resizer(this->flags_handler_,n_basis,n_points);

        // apply the functor resizer to each pair ValueType/Container in this->values_
        boost::fusion::for_each(this->values_,resizer);

        this->set_initialized(true);
    }

};


template<int dim, int codim, int range, int rank>
class FuncValuesCache : public ValuesCache<dim,codim,range,rank,ValueVector>
{
    using parent_t = ValuesCache<dim,codim,range,rank,ValueVector>;

private:
    struct FuncValuesCacheResizer
    {
        FuncValuesCacheResizer(
            FunctionFlags &flags_handler,
            const int n_points)
            :
            flags_handler_ {flags_handler},
        n_points_ {n_points}
        {
            Assert(n_points >= 0, ExcLowerRange(n_points,1));
        }

        template <class pair_ValueType_ValueContainer>
        void operator()(pair_ValueType_ValueContainer &x) const
        {
            using ValueType = typename pair_ValueType_ValueContainer::first_type;
            auto &value = x.second;

            if (this->flags_handler_.template fill<ValueType>())
            {
                if (value.get_num_points() != n_points_)
                {
                    value.resize(n_points_);
                }
                value.zero();
            }
            else
            {
                value.clear();
                this->flags_handler_.template set_filled<ValueType>(false);
            }
        }

    private:
        FunctionFlags &flags_handler_;

        int n_points_;
    };



public:
    using Func = Function<dim,codim,range,rank>;

    using Point = typename Func::Point;

    /**
     * Allocate space for the fucntion values and derivatives
     * of the element at quadrature points
     * as specify by the flag
     */
    void resize(const FunctionFlags &flags_handler,
                const Size n_points)
    {
        this->flags_handler_ = flags_handler;

        FuncValuesCacheResizer resizer(this->flags_handler_,n_points);
        // apply the functor resizer to each pair ValueType/Container in this->values_
        boost::fusion::for_each(this->values_,resizer);

        this->set_initialized(true);
    }

    ValueVector<Point> points_;

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
        return std::get<topology_dim>(values_)[topology_id];
    }

    template <int topology_dim>
    const Cache &
    get_value_cache(const int topology_id) const
    {
        const auto &cache = std::get<topology_dim>(values_)[topology_id];
        Assert(cache.is_filled() == true, ExcCacheNotFilled());

        return cache;
    }

    CacheList<Cache, Cache::get_dim()> values_;

};


IGA_NAMESPACE_CLOSE



#endif // #ifndef VALUES_CACHE_H_

