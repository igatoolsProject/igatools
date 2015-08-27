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

#ifndef __FUNCTION_ELEMENT_H
#define __FUNCTION_ELEMENT_H

#include <igatools/geometry/grid_element.h>
#include <igatools/base/value_types.h>
#include <igatools/basis_functions/values_cache.h>

IGA_NAMESPACE_OPEN

namespace function_element
{
enum class Flags
{
    /** Fill nothing */
    none           =    0,

    /** Quadrature points on the element */
    value          =    1L << 1,

    /** Quadrature weigths on the element */
    gradient       =    1L << 2
};

template <int,int,int,int> class Function;
template <int,int> class PhysicalDomain;
template <int,int> class PhysicalDomainElement;

/**
 *
 * @ingroup serializable
 */
template<int dim, int codim, int range, int rank, class ContainerType_>
class FunctionElementBase
{
private:
    using self_t = FunctionElementBase<dim, codim, range, rank, ContainerType_>;

public:

    using ContainerType = ContainerType_;
    using Func = ContainerType_;
    using MapFunc = typename Func::MapFunc;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Gradient = typename Func::Gradient;
    using Hessian  = typename Func::Hessian;
    using Div      = typename Func::Div;

    using Flags = function_element::Flags;

    using Grid = CartesianGrid<dim>;
    //using IndexType = typename Grid::IndexType;
    //using List = typename Grid::List;
    using ListIt = typename Grid::ListIt;

private:
    template <int order>
    using Derivative = typename Func::template Derivative<order>;

protected:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    FunctionElementBase() = default;

public:
    /**
     * Construct an accessor pointing to the element with
     * flat index @p elem_index of the CartesianGrid @p grid.
     */
    FunctionElementBase(const std::shared_ptr<ContainerType_> func,
                        const ListIt &index,
                        const PropId &prop = ElementProperties::active);

    /**
     * Copy constructor.
     * It can be used with different copy policies
     * (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a
     * classic copy constructor)
     * uses the deep copy.
     */
    FunctionElementBase(const self_t &elem,
                        const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    FunctionElement(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~FunctionElement() = default;
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
    void deep_copy_from(const self_t &element)
    {
        Assert(false,ExcNotImplemented());
    }

    /**
     * Performs a shallow copy of the input @p element. The current object will contain a pointer to the
     * local cache used by the input @p element.
     */
    void shallow_copy_from(const self_t &element)
    {
        Assert(false,ExcNotImplemented());
    }

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



    template<class ValueType, int sdim>
    auto
    get_values(const int s_id) const
    {
        Assert(all_sub_elems_cache_ != nullptr,ExcNullPtr());
        const auto &cache =
            all_sub_elems_cache_->template get_sub_elem_cache<sdim>(s_id);
        return cache.template get_data<ValueType>();
    }


#if 0
    /**
     * @name Methods for the for the evaluations of Functions's derivatives
     *  without the use of the cache.
     */
    ///@{
    /**
     * Returns a ValueTable with the values specified by the template parameter
     * <tt>ValueType</tt>
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_cache()/fill_cache().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    template <class ValueType>
    decltype(auto) evaluate_at_points(const Quadrature<dim> &points)
    {
        func_->reset(ValueType::flag,points);
        const auto topology = Topology<dim>();
        func_->init_cache(*this,topology);
        func_->fill_cache(*this,topology,0);

        return this->template get_values<ValueType,dim>(0);
    }
    ///@}


    /**
     * Returns the flags that are valid to be used with this class.
     *
     * @note The valid flags are defined to be the ones that can be inferred from the ValueType(s)
     * used as key of the boost::fusion::map in CType.
     */
    static ValueFlags get_valid_flags();
#endif

    /**
     * Returns the <tt>topology_dim</tt> dimensional topology_id-th sub-element measure
     * multiplied by the weights of the quadrature.
     */
    template <int topology_dim>
    ValueVector<Real> get_w_measures(const int topology_id) const
    {
        ValueVector<Real> w_meas;
        Assert(false,ExcNotImplemented());

        return w_meas;
    }



    typename List::iterator &operator++()
    {
        return (++(*phys_domain_elem_));
    }

private:

    using CType = boost::fusion::map<
                  boost::fusion::pair<     _Value,DataWithFlagStatus<ValueVector<Value>>>,
                  boost::fusion::pair<  _Gradient,DataWithFlagStatus<ValueVector<Derivative<1>>>>,
                  boost::fusion::pair<   _Hessian,DataWithFlagStatus<ValueVector<Derivative<2>>>>,
                  boost::fusion::pair<_Divergence,DataWithFlagStatus<ValueVector<Div>>>,
                  boost::fusion::pair<     _Point,DataWithFlagStatus<ValueVector<Point>>>
                  >;



    using Cache = FuncValuesCache<dim,CType>;




private:
    using CacheType = AllSubElementsCache<Cache>;

    //TODO (martinelli, Aug 13, 2015): this function should not be public.
    std::shared_ptr<CacheType>
    &get_cache()
    {
        Assert(all_sub_elems_cache_ != nullptr,ExcNullPtr());
        return all_sub_elems_cache_;
    }

private:
    //using PhysDomain = typename Func::PhysDomain;
    using PhysDomainElem = typename Func::PhysDomain::ElementAccessor;

    std::shared_ptr<ContainerType_> func_;

    std::shared_ptr<PhysDomainElem> phys_domain_elem_;

    std::shared_ptr<AllSubElementsCache<Cache>> all_sub_elems_cache_;

    template <class Accessor> friend class GridIteratorBase;
    friend class Function<dim, codim, range, rank>;

public:
    const PhysicalDomainElement &get_domain_element() const;

    void print_info(LogStream &out) const;

    void print_cache_info(LogStream &out) const;


    /**
     * @name Comparison operators.
     *
     * @brief The comparison operators compares the <em>position</em> of the element in the grid.
     *
     * @warning To be comparable, two Function objects must be defined using the same Function
     * (and therefore on the same grid),
     * otherwise an assertion will be raised (in Debug mode).
     */
    ///@{
    /** Returns TRUE if the two elements have the same index on the grid. */
    bool operator==(const self_t &a) const;


    /** Returns TRUE if the two elements have different indices on the grid. */
    bool operator!=(const self_t &a) const;

    /**
     * Returns TRUE if the the index of the element on the left of the operator <tt> < </tt>
     * is smaller than the the index of the element on the right.
     * */
    bool operator<(const self_t &a) const;

    /**
     * Returns TRUE if the the index of the element on the left of the operator <tt> < </tt>
     * is bigger than the the index of the element on the right.
     * */
    bool operator>(const self_t &a) const;
    ///@}

//
//    /** Returns the index of the element. */
//    IndexType get_index() const;



private:
#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif // SERIALIZATION
};


template <int dim, int codim, int range, int rank>
class ConstFunctionElement
    : public FunctionElementBase<dim, codim, range, rank,
      const Function<dim,codim,range,rank>>
{
    using FunctionElementBase<dim, codim, range, rank,
          const Function<dim,codim,range,rank>>::FunctionElementBase;
};


//template <int dim, int codim, int range, int rank>
//class FunctionElement
//    : public FunctionElementBase<dim, CartesianFunction<dim>>
//{
//    using FunctionElementBase<dim, CartesianFunction<dim>>::FunctionElementBase;
//public:
//    void add_property(const PropId &prop)
//    {
//        this->grid_->elem_properties_[prop].insert(this->get_index());
//    }
//};


IGA_NAMESPACE_CLOSE

#endif
