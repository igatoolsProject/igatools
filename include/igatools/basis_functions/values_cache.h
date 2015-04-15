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
//#include <igatools/utils/static_multi_array.h>
//#include <igatools/utils/cartesian_product_indexer.h>

//#include <igatools/basis_functions/spline_space.h>

//#include <igatools/basis_functions/space_element_base.h>



IGA_NAMESPACE_OPEN


template<int dim_, int codim, int range, int rank, template<class ContainedType> class ContainerType>
class ValuesCache : public CacheStatus
{
    using Func = Function<dim_,codim,range,rank>;

    using Value = typename Func::Value;

    template <int der_order>
    using Derivative = typename Func::template Derivative<der_order>;

    using Div = typename Func::Div;

public:

    static constexpr int dim = dim_;





    void print_info(LogStream &out) const
    {
        out.begin_item("Fill flags:");
        flags_handler_.print_info(out);
        out.end_item();

        if (flags_handler_.template filled<_Value>())
        {
            out.begin_item("Values:");
            get_der<_Value>().print_info(out);
            out.end_item();
        }

        if (flags_handler_.template filled<_Gradient>())
        {
            out.begin_item("Gradients:");
            get_der<_Gradient>().print_info(out);
            out.end_item();
        }

        if (flags_handler_.template filled<_Hessian>())
        {
            out.begin_item("Hessians:");
            get_der<_Hessian>().print_info(out);
            out.end_item();
        }

        if (flags_handler_.template filled<_Divergence>())
        {
            out.begin_item("Divergences:");
            get_der<_Divergence>().print_info(out);
            out.end_item();
        }
    }


    FunctionFlags flags_handler_;


    using map_TuplePosition_ContainerType = boost::mpl::map<
                                            boost::mpl::pair<TuplePosition_from_ValueType<     _Value>,ContainerType<Value> >,
                                            boost::mpl::pair<TuplePosition_from_ValueType<  _Gradient>,ContainerType<Derivative<1>> >,
                                            boost::mpl::pair<TuplePosition_from_ValueType<   _Hessian>,ContainerType<Derivative<2>> >,
                                            boost::mpl::pair<TuplePosition_from_ValueType<_Divergence>,ContainerType<Div>>
                                            >;
    using map_TP_CT = map_TuplePosition_ContainerType;

    template <int tuple_position>
    using ContType_from_TuplePos = typename boost::mpl::at<map_TP_CT,boost::mpl::int_<tuple_position>>::type;

    std::tuple<
    ContType_from_TuplePos<0>,
                           ContType_from_TuplePos<1>,
                           ContType_from_TuplePos<2>,
                           ContType_from_TuplePos<3>> values_;

    template<class ValueType>
    auto &get_der()
    {
        return std::get<TuplePosition_from_ValueType<ValueType>::value>(values_);
    }

    template<class ValueType>
    const auto &get_der() const
    {
        //TODO (martinelli, Apr 03,2015): uncomment this assertion
//        Assert(flags_handler_.filled<ValueType>(),
//               ExcMessage("The cache for " + ValueType::name + " is not filled."));

        return std::get<TuplePosition_from_ValueType<ValueType>::value>(values_);
    }



    template<class ValueType>
    void clear_der()
    {
        auto &value = std::get<TuplePosition_from_ValueType<ValueType>::value>(values_);
        value.clear();
    }


};




template<int dim, int codim, int range, int rank>
class BasisValuesCache : public ValuesCache<dim,codim,range,rank,ValueTable>
{
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

        Assert(n_points >= 0, ExcLowerRange(n_points,1));
        Assert(n_basis > 0, ExcLowerRange(n_basis,1));

        if (this->flags_handler_.template fill<_Value>())
            resize_der<_Value>(n_basis,n_points);
        if (this->flags_handler_.template fill<_Gradient>())
            resize_der<_Gradient>(n_basis,n_points);
        if (this->flags_handler_.template fill<_Hessian>())
            resize_der<_Hessian>(n_basis,n_points);

        if (this->flags_handler_.template fill<_Divergence>())
        {
            Assert(this->flags_handler_.template fill<_Gradient>(),
                   ExcMessage("Divergence requires gradient to be filled."));
            resize_der<_Divergence>(n_basis,n_points);
        }

        this->set_initialized(true);
    }

    template<class ValueType>
    void resize_der(const int n_basis, const int n_points)
    {
        auto &value = this->template get_der<ValueType>();
        if (value.get_num_points() != n_points ||
            value.get_num_functions() != n_basis)
        {
            value.resize(n_basis, n_points);
            value.zero();
        }
    }

};


template<int dim, int codim, int range, int rank>
class FuncValuesCache : public ValuesCache<dim,codim,range,rank,ValueVector>
{
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

        Assert(n_points >= 0, ExcLowerRange(n_points,1));

        if (this->flags_handler_.fill_points())
            points_.resize(n_points);

        if (this->flags_handler_.template fill<_Value>())
            resize_der<_Value>(n_points);
        if (this->flags_handler_.template fill<_Gradient>())
            resize_der<_Gradient>(n_points);
        if (this->flags_handler_.template fill<_Hessian>())
            resize_der<_Hessian>(n_points);

        if (this->flags_handler_.template fill<_Divergence>())
        {
            Assert(this->flags_handler_.template fill<_Gradient>(),
                   ExcMessage("Divergence requires gradient to be filled."));

            Assert(false,ExcNotImplemented());
            AssertThrow(false,ExcNotImplemented());
            resize_der<_Divergence>(n_points);
        }

        this->set_initialized(true);
    }

    template<class ValueType>
    void resize_der(const int n_points)
    {
        auto &value = this->template get_der<ValueType>();
        if (value.get_num_points() != n_points)
        {
            value.resize(n_points);
            value.zero();
        }
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
    get_value_cache(const int j)
    {
        return std::get<topology_dim>(values_)[j];
    }

    template <int topology_dim>
    const Cache &
    get_value_cache(const int j) const
    {
        return std::get<topology_dim>(values_)[j];
    }

//        CacheList<ValuesCache, dim> values_;

    CacheList<Cache, Cache::dim> values_;

};


IGA_NAMESPACE_CLOSE



#endif // #ifndef VALUES_CACHE_H_

