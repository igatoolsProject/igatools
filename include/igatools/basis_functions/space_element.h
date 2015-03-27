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


#ifndef SPACE_ELEMENT_H_
#define SPACE_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/base/value_types.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/flags_handler.h>

#include <igatools/base/function.h>

#include <igatools/base/quadrature.h>

#include <igatools/utils/value_vector.h>
#include <igatools/utils/value_table.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/cartesian_product_indexer.h>

#include <igatools/basis_functions/spline_space.h>

#include <igatools/basis_functions/space_element_base.h>


#include <boost/mpl/map.hpp>
#include <boost/mpl/int.hpp>
//#include <boost/mpl/order.hpp>

IGA_NAMESPACE_OPEN




template<int dim,int codim,int range,int rank>
class SpaceElement : public SpaceElementBase<dim>
{
protected:
    using base_t =  SpaceElementBase<dim>;
private:
    using self_t = SpaceElement<dim,codim,range,rank>;

public:

    using Func = Function<dim,codim,range,rank>;

    using RefPoint = typename Func::RefPoint;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    template <int order>
    using Derivative = typename Func::template Derivative<order>;
    using Div = typename Func::Div;

    static const int space_dim = Func::space_dim;

    using Topology = typename base_t::Topology;

    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentContainer = typename SplineSpace<dim,range,rank>::template ComponentContainer<T>;
    using TensorSizeTable = typename SplineSpace<dim,range,rank>::TensorSizeTable;
    ///@}


    /** @name Constructors */
    ///@{
    using SpaceElementBase<dim>::SpaceElementBase;

    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a classic copy constructor)
     * uses the deep copy.
     */
    SpaceElement(const self_t &elem,
                 const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    SpaceElement(self_t &&elem) = default;

    /**
     * Destructor.
     */
    virtual ~SpaceElement() = default;
    ///@}

    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator. Performs a <b>shallow copy</b> of the input @p element.
     *
     * @note Internally it uses the function shallow_copy_from().
     */
    self_t &operator=(const self_t &element);

    /**
     * Move assignment operator.
     */
    self_t &operator=(self_t &&elem) = default;
    ///@}

    /**
     * @name Functions for performing different kind of copy.
     */
    ///@{
    /**
     * Performs a deep copy of the input @p element,
     * i.e. a new local cache is built using the copy constructor on the local cache of @p element.
     *
     * @note In DEBUG mode, an assertion will be raised if the input local cache is not allocated.
     */
    void deep_copy_from(const self_t &element);

    /**
     * Performs a shallow copy of the input @p element. The current object will contain a pointer to the
     * local cache used by the input @p element.
     */
    void shallow_copy_from(const self_t &element);
    ///@}


    template <class ValueType, int k = dim>
    auto
    get_basis(const int j, const std::string &dofs_property = DofProperties::active) const
    {
        Assert(local_cache_ != nullptr, ExcNullPtr());
        const auto &cache = local_cache_->template get_value_cache<k>(j);
        Assert(cache.is_filled() == true, ExcCacheNotFilled());
        const auto values_all_elem_dofs = cache.template get_der<ValueType>();

        //--------------------------------------------------------------------------------------
        // filtering the values that correspond to the dofs with the given property --- begin
        vector<Index> dofs_global;
        vector<Index> dofs_local_to_patch;
        vector<Index> dofs_local_to_elem;

        this->space_->get_element_dofs(
            this->as_cartesian_grid_element_accessor(),
            dofs_global,
            dofs_local_to_patch,
            dofs_local_to_elem,
            dofs_property);

        const auto n_filtered_dofs = dofs_local_to_elem.size();
        const auto n_pts = values_all_elem_dofs.get_num_points();

        decltype(values_all_elem_dofs) values_filtered_elem_dofs(n_filtered_dofs,n_pts);

        int fn = 0;
        for (const auto loc_dof : dofs_local_to_elem)
        {
            const auto values_all_elem_dofs_fn = values_all_elem_dofs.get_function_view(loc_dof);

            const auto values_filtered_elem_dofs_fn = values_filtered_elem_dofs.get_function_view(fn);

            std::copy(values_all_elem_dofs_fn.begin(),
                      values_all_elem_dofs_fn.end(),
                      values_filtered_elem_dofs_fn.begin());

            ++fn;
        }
        // filtering the values that correspond to the dofs with the given property --- end
        //--------------------------------------------------------------------------------------

        return values_filtered_elem_dofs;
//        return this->template get_values<ValueType::order,k>(j,dofs_property);
    }

#if 0
    template<int order = 0, int k = dim>
    auto
    get_values(const int j, const std::string &dofs_property = DofProperties::active) const
    {
        Assert(local_cache_ != nullptr, ExcNullPtr());
        const auto &cache = local_cache_->template get_value_cache<k>(j);
        Assert(cache.is_filled() == true, ExcCacheNotFilled());
        const auto values_all_elem_dofs = cache.template get_der<order>();

        //--------------------------------------------------------------------------------------
        // filtering the values that correspond to the dofs with the given property --- begin
        vector<Index> dofs_global;
        vector<Index> dofs_local_to_patch;
        vector<Index> dofs_local_to_elem;

        this->space_->get_element_dofs(
            this->as_cartesian_grid_element_accessor(),
            dofs_global,
            dofs_local_to_patch,
            dofs_local_to_elem,
            dofs_property);

        const auto n_filtered_dofs = dofs_local_to_elem.size();
        const auto n_pts = values_all_elem_dofs.get_num_points();

        decltype(values_all_elem_dofs) values_filtered_elem_dofs(n_filtered_dofs,n_pts);

        int fn = 0;
        for (const auto loc_dof : dofs_local_to_elem)
        {
            const auto values_all_elem_dofs_fn = values_all_elem_dofs.get_function_view(loc_dof);

            const auto values_filtered_elem_dofs_fn = values_filtered_elem_dofs.get_function_view(fn);

            std::copy(values_all_elem_dofs_fn.begin(),
                      values_all_elem_dofs_fn.end(),
                      values_filtered_elem_dofs_fn.begin());

            ++fn;
        }
        // filtering the values that correspond to the dofs with the given property --- end
        //--------------------------------------------------------------------------------------

        return values_filtered_elem_dofs;
    }
#endif

    auto
    get_element_values(const std::string &dofs_property = DofProperties::active) const
    {
        return this->template get_basis<_Value,dim>(0,dofs_property);
//        return this->template get_values<0,dim>(0,dofs_property);
    }

    template <class ValueType, int k = dim>
    auto
    linear_combination(const vector<Real> &loc_coefs,
                       const int id,
                       const std::string &dofs_property) const
    {
        const auto &basis_values =
            this->template get_basis<ValueType, k>(id,dofs_property);
        return basis_values.evaluate_linear_combination(loc_coefs) ;
    }


    template<int k = dim>
    ValueTable<Div> get_divergences(const int id,
                                    const std::string &dofs_property) const
    {
        /*
        Assert(local_cache_ != nullptr, ExcNullPtr());
        const auto &cache = local_cache_->template get_value_cache<k>(id);
        Assert(cache.is_filled() == true, ExcCacheNotFilled());

        Assert(cache.flags_handler_.gradients_filled() == true, ExcCacheNotFilled());
        //*/
        const auto &basis_gradients =
            this->template get_basis<_Gradient,k>(id,dofs_property);

        const int n_basis = basis_gradients.get_num_functions();
        const int n_pts   = basis_gradients.get_num_points();

        ValueTable<Div> divergences(n_basis,n_pts);
        /*
        std::transform(basis_gradients.cbegin(),
                       basis_gradients.cend(),
                       divergences.begin(),
                       [](const auto &grad){ return trace(grad);});
                       //*/

        auto div_it = divergences.begin();
        for (const auto &grad : basis_gradients)
        {
            *div_it = trace(grad);
            ++div_it;
        }
//*/
        return divergences;
    }



    ValueTable<Div> get_element_divergences(const std::string &dofs_property) const
    {
        return get_divergences<dim>(0,dofs_property);
    }




    void print_info(LogStream &out) const;

    void print_cache_info(LogStream &out) const;

protected:


    /**
     * Base class for the cache of the element values and
     * for the cache of the face values.
     */
    class ValuesCache : public CacheStatus
    {
    public:
        /**
         * Allocate space for the values and derivatives
         * of the element basis functions at quadrature points
         * as specify by the flag
         */
        void resize(const FunctionFlags &flags_handler,
                    const Size total_n_points,
                    const Size n_basis);



    public:
        void print_info(LogStream &out) const;

        FunctionFlags flags_handler_;

        using map_ValueType_TuplePosition = boost::mpl::map<
                                            boost::mpl::pair<     _Value,boost::mpl::int_<0> >,
                                            boost::mpl::pair<  _Gradient,boost::mpl::int_<1> >,
                                            boost::mpl::pair<   _Hessian,boost::mpl::int_<2> >,
                                            boost::mpl::pair<_Divergence,boost::mpl::int_<3> >
                                            >;
        using map_VT_TP = map_ValueType_TuplePosition;

        template <class ValueType>
        using TuplePos_from_ValueType = typename boost::mpl::at<map_VT_TP,ValueType>::type;

        using map_TuplePosition_ContainerType = boost::mpl::map<
                                                boost::mpl::pair<TuplePos_from_ValueType<     _Value>,ValueTable<Value> >,
                                                boost::mpl::pair<TuplePos_from_ValueType<  _Gradient>,ValueTable<Derivative<1>> >,
                                                boost::mpl::pair<TuplePos_from_ValueType<   _Hessian>,ValueTable<Derivative<2>> >,
                                                boost::mpl::pair<TuplePos_from_ValueType<_Divergence>,ValueVector<Derivative<1>>>
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
            return std::get<TuplePos_from_ValueType<ValueType>::value>(values_);
        }

        template<class ValueType>
        const auto &get_der() const
        {
            const auto k = TuplePos_from_ValueType<ValueType>::value;

            Assert(flags_handler_.filled<ValueType>(),
                   ExcMessage(ValueType::name + "s cache is not filled."));

            return std::get<k>(values_);
        }


        template<class ValueType>
        void resize_der(const int n_basis, const int n_points)
        {
            auto &value = std::get< TuplePos_from_ValueType<ValueType>::value >(values_);
            if (value.get_num_points() != n_points ||
                value.get_num_functions() != n_basis)
            {
                value.resize(n_basis, n_points);
                value.zero();
            }
        }

        template<class ValueType>
        void clear_der()
        {
            auto &value = std::get<ValueType::order>(values_);
            value.clear();
        }

    };


    class LocalCache
    {
    public:
        LocalCache() = default;

        LocalCache(const LocalCache &in) = default;

        LocalCache(LocalCache &&in) = default;

        ~LocalCache() = default;


        LocalCache &operator=(const LocalCache &in) = delete;

        LocalCache &operator=(LocalCache &&in) = delete;

        void print_info(LogStream &out) const;

        template <int topology_dim>
        ValuesCache &
        get_value_cache(const int j)
        {
            return std::get<topology_dim>(values_)[j];
        }

        template <int topology_dim>
        const ValuesCache &
        get_value_cache(const int j) const
        {
            return std::get<topology_dim>(values_)[j];
        }

        CacheList<ValuesCache, dim> values_;
    };

    /** The local (element and face) cache. */
    std::shared_ptr<LocalCache> local_cache_;

public:
    // TODO (pauletti, Mar 17, 2015): this cannot be public, if needed it means wrong desing
    std::shared_ptr<LocalCache> &get_local_cache()
    {
        return local_cache_;
    }

protected:
    /**
     * Performs a copy of the input @p element.
     * The type of copy (deep or shallow) is specified by the input parameter @p copy_policy.
     */
    void copy_from(const self_t &element,
                   const CopyPolicy &copy_policy);




};


IGA_NAMESPACE_CLOSE



#endif // #ifndef SPACE_ELEMENT_H_

