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

#ifndef FUNCTION_ELEMENT_H
#define FUNCTION_ELEMENT_H

#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/base/value_types.h>
#include <igatools/basis_functions/values_cache.h>

IGA_NAMESPACE_OPEN

template <int,int,int,int> class Function;

template <int,int> class Mapping;
template <int,int> class MappingElement;

/**
 *
 * @ingroup serializable
 */
template<int dim, int codim, int range = 1, int rank = 1>
class FunctionElement
{
private:
	using self_t = FunctionElement<dim,codim,range,rank>;

public:
    using Func = Function<dim, codim, range, rank>;
    using MapFunc = typename Func::MapFunc;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Gradient = typename Func::Gradient;
    using Hessian  = typename Func::Hessian;
    using Div      = typename Func::Div;
//    using ContainerType = const CartesianGrid<dim>;
    using ContainerType = const Func;

private:
    template <int order>
    using Derivative = typename Func::template Derivative<order>;

public:

    /** @name Constructors */
    ///@{
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    FunctionElement() = default;


    /**
     * Construct an accessor pointing to the element with
     * flat index @p elem_index of the CartesianGrid @p grid.
     */
    FunctionElement(const std::shared_ptr<const Func> func,
                    const Index elem_index);

    /**
     * Copy constructor.
     * It can be used with different copy policies
     * (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a
     * classic copy constructor)
     * uses the deep copy.
     */
    FunctionElement(const FunctionElement<dim,codim,range,rank> &elem,
                    const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    FunctionElement(FunctionElement<dim,codim,range,rank> &&elem) = default;

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
    void deep_copy_from(const FunctionElement<dim,codim,range,rank> &element)
    {
        Assert(false,ExcNotImplemented());
    }

    /**
     * Performs a shallow copy of the input @p element. The current object will contain a pointer to the
     * local cache used by the input @p element.
     */
    void shallow_copy_from(const FunctionElement<dim,codim,range,rank> &element)
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
    FunctionElement<dim,codim,range,rank> &operator=(const FunctionElement<dim,codim,range,rank> &element);

    /**
     * Move assignment operator.
     */
    FunctionElement<dim,codim,range,rank> &operator=(FunctionElement<dim,codim,range,rank> &&elem) = default;
    ///@}



    template<class ValueType, int k>
    auto
    get_values(const int j) const
    {
        Assert(all_sub_elems_cache_ != nullptr,ExcNullPtr());
        const auto &cache = all_sub_elems_cache_->template get_sub_elem_cache<k>(j);
        return cache.template get_data<ValueType>();
    }



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

private:

    using CType = boost::fusion::map<
                  boost::fusion::pair<     _Value,DataWithFlagStatus<ValueVector<Value>>>,
                  boost::fusion::pair<  _Gradient,DataWithFlagStatus<ValueVector<Derivative<1>>>>,
                  boost::fusion::pair<   _Hessian,DataWithFlagStatus<ValueVector<Derivative<2>>>>,
                  boost::fusion::pair<_Divergence,DataWithFlagStatus<ValueVector<Div>>>,
                  boost::fusion::pair<     _Point,DataWithFlagStatus<ValueVector<Point>>>
                  >;



    using Cache = FuncValuesCache<dim,CType>;

    std::shared_ptr<AllSubElementsCache<Cache>> all_sub_elems_cache_;


public:
    using CacheType = AllSubElementsCache<Cache>;

    //TODO (martinelli, Aug 13, 2015): this function should not be public.
    std::shared_ptr<CacheType>
    &get_cache()
	{
        Assert(all_sub_elems_cache_ != nullptr,ExcNullPtr());
    	return all_sub_elems_cache_;
	}

private:


    std::shared_ptr<Func> func_;

    using GridElem = CartesianGridElement<dim>;
    std::shared_ptr<GridElem> grid_elem_;


    using PhysDomain = Mapping<dim,codim>;
    using PhysDomainElem = MappingElement<dim,codim>;
    std::shared_ptr<PhysDomainElem> phys_domain_elem_;

    template <class Accessor> friend class CartesianGridIteratorBase;
    friend class Function<dim, codim, range, rank>;

    /**
     * Creates a new object performing a deep copy of the current object using the FunctionElement
     * copy constructor.
     */
    std::shared_ptr<FunctionElement<dim,codim,range,rank> > clone() const;



public:
    const GridElem & get_grid_element() const;

    void print_info(LogStream &out) const;

    void print_cache_info(LogStream &out) const;

    std::shared_ptr<const CartesianGrid<dim>> get_grid() const;

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

    /**
     * Sets the index of the element using the flatten representation.
     * @note This function also updates the index for the tensor representation.
     * @warning This may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    void move_to(const Index flat_index) ;


    /** @name Functions related to the indices of the element in the cartesian grid. */
    ///@{
    /** Returns the index of the element in its flatten representation. */
    Index get_flat_index() const;

    /** Returns the index of the element in its tensor representation. */
    TensorIndex<dim> get_tensor_index() const;
    ///@}



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


#if 0
template<int dim, int codim, int range = 1, int rank = 1>
class IgFunctionElement
    : public FunctionElement<dim,codim,range,rank>
{

};
#endif

IGA_NAMESPACE_CLOSE

#endif
