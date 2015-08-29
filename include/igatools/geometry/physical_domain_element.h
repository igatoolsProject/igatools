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

#ifndef __PHYSICAL_DOMAIN_ELEMENT_H_
#define __PHYSICAL_DOMAIN_ELEMENT_H_

#include <igatools/geometry/grid_element.h>
#include <igatools/geometry/physical_domain.h>

IGA_NAMESPACE_OPEN

template <int,int,int,int> class Function;
template <int,int,int,int, class> class FunctionElementBase;
template <int,int,int,int> class FunctionElement;
//template <int,int> class PhysicalDomain;

/**
 *
 * @ingroup elements
 */
template<int dim_, int codim_, class ContainerType_>
class PhysicalDomainElementBase
{
private:
    using self_t  = PhysicalDomainElementBase<dim_, codim_, ContainerType_>;

public:
    using ContainerType = ContainerType_;
    using GridElem = typename ContainerType_::GridType::ConstElementAccessor;
    using FuncElem = typename ContainerType_::FuncType::ConstElementAccessor;
    using ListIt = typename ContainerType_::ListIt;

    using Flags = physical_domain_element::Flags;

    // using FuncType = typename ContainerType_::FuncType;

//    static const int dim = dim_;
//    static const int codim = codim_;
//    static const int space_dim = dim_+codim_;

    //using ListIt = typename FuncType::ListIt;
    //using FuncElem = typename FuncType::ElementAccessor;


    /** @name Constructors */
    ///@{
protected:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    PhysicalDomainElementBase() = default;

public:
    /**
     * Construct an accessor pointing to the element with
     * flat index @p elem_index of the Function @p func.
     */
    PhysicalDomainElementBase(const std::shared_ptr<ContainerType_> phys_dom,
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
    PhysicalDomainElementBase(const self_t &elem,
                              const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    PhysicalDomainElementBase(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~PhysicalDomainElementBase() = default;
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

    /**
     * @name Comparison operators
     * @note In order to be meaningful, the comparison must be performed on elements defined on
     * the <b>same grid</b>
     * (in the sense that the pointer to the grid held by the element must point to the same
     * grid object).
     */
    ///@{
    /**
     * True if the elements have the same index.
     *  @note In debug mode, it is also check they both refer to
     *  the same cartesian grid. No check is done on the cache.
     */
    bool operator==(const self_t &elem) const;

    /**
     * True if the elements have different index.
     *  @note In debug mode, it is also check they both refer to
     *  the same cartesian grid. No check is done on the cache.
     */
    bool operator!=(const self_t &elem) const;

    /**
     * True if the flat-index of the element on the left is smaller than
     * the flat-index of the element on the right.
     *  @note In debug mode, it is also check they both refer to
     *  the same cartesian grid. No check is done on the cache.
     */
    bool operator<(const self_t &elem) const;

    /**
     * True if the flat-index of the element on the left is bigger than
     * the flat-index of the element on the right.
     *  @note In debug mode, it is also check they both refer to
     *  the same cartesian grid. No check is done on the cache.
     */
    bool operator>(const self_t &elem) const;
    ///@}

#if 0
    template<int order>
    using InvDerivative = typename FuncType::template InvDerivative<order>;

    template <int order>
    using Derivative = typename FuncType::template Derivative<order>;
#endif

public:
    ListIt &operator++()
    {
        return (++(*grid_elem_));
    }


    std::shared_ptr<GridElem> get_grid_element() const
    {
        return grid_elem_;
    }

public:
    template<int sdim>
    ValueVector<Real> const &get_measures(const int s_id) const
    {
        return get_values_from_cache<_Measure,sdim>(s_id);
    }

    template<int sdim>
    ValueVector<Real> const &get_w_measures(const int s_id) const
    {
        return get_values_from_cache<_W_Measure,sdim>(s_id);
    }
#if 0
    const ValueVector<Points<space_dim> > &get_external_normals() const;

    using MetricTensor =
        Tensor<dim, 1, tensor::covariant, Tensor<dim, 1, tensor::contravariant, Tdouble> >;

    ValueVector<MetricTensor> compute_inv_first_fundamental_form() const;

    ValueVector<MetricTensor> compute_second_fundamental_form() const;

    ValueVector< Derivative<1> > get_D_external_normals() const;

    const ValueVector<SafeSTLVector<Real> > &get_principal_curvatures() const;


    template<int sub_dim>
    const ValueVector<Points<space_dim> > &
    get_boundary_normals(const int s_id) const
    {
#if 0
        Assert(dim==sub_dim+1, ExcNotImplemented());
        ValueVector<Points<space_dim>> res;
        const auto &DF_inv = get_values_from_cache<_InvGradient, sub_dim>(s_id);
        const auto n_hat  = this->get_grid()->template get_boundary_normals<sub_dim>(s_id)[0];

        const auto n_points = DF_inv.get_num_points();
        res.resize(n_points);
        for (int pt = 0; pt < n_points; ++pt)
        {
            const auto DF_inv_t = co_tensor(transpose(DF_inv[pt]));
            res[pt] = action(DF_inv_t, n_hat);
            res[pt] /= res[pt].norm();
        }
        return res;
#endif
        return get_values_from_cache<_BoundaryNormal,sub_dim>(s_id);
    }
#endif

private:
    template <class ValueType, int topology_dim = dim_>
    auto &get_values_from_cache(const int topology_id = 0) const
    {
        Assert(local_cache_ != nullptr,ExcNullPtr());
        const auto &cache = local_cache_->template
                            get_sub_elem_cache<topology_dim>(topology_id);
        return cache.template get_data<ValueType>();
    }
private:

    using CType = boost::fusion::map<
                  boost::fusion::pair<       _Measure,DataWithFlagStatus<ValueVector<Real>>>,
                  boost::fusion::pair<     _W_Measure,DataWithFlagStatus<ValueVector<Real>>>>;
//                ,
//                  boost::fusion::pair<   _InvGradient,DataWithFlagStatus<ValueVector<InvDerivative<1>>>>,
//                  boost::fusion::pair<    _InvHessian,DataWithFlagStatus<ValueVector<InvDerivative<2>>>>,
//                  boost::fusion::pair<_BoundaryNormal,DataWithFlagStatus<ValueVector<Points<space_dim>>>>,
//                  boost::fusion::pair<   _OuterNormal,DataWithFlagStatus<ValueVector<Points<space_dim>>>>,
//                  boost::fusion::pair<     _Curvature,DataWithFlagStatus<ValueVector<SafeSTLVector<Real>>>>
//                  >;


//    /**
//     * Returns the flags that are valid to be used with this class.
//     *
//     * @note The valid flags are defined to be the ones that can be inferred from the ValueType(s)
//     * used as key of the boost::fusion::map in CType.
//     */
//    static ValueFlags get_valid_flags()
//    {
//        return cacheutils::get_valid_flags_from_cache_type(CType());
//    }

    using Cache = FuncValuesCache<dim_,CType>;


public:
    using CacheType = AllSubElementsCache<Cache>;

private:


    std::shared_ptr<ContainerType_> phys_dom_;

    std::shared_ptr<GridElem> grid_elem_;

    std::shared_ptr<FuncElem> func_elem_;

    std::shared_ptr<CacheType> local_cache_;

    template <class Accessor> friend class GridIteratorBase;
    friend class PhysicalDomain<dim_, codim_>;

//    /**
//     * Creates a new object performing a deep copy of the current object using the PhysicalDomainElement
//     * copy constructor.
//     */
//    std::shared_ptr<PhysicalDomainElement<dim_,codim_> > clone() const;

};


template <int dim, int codim>
class ConstPhysicalDomainElement
    : public PhysicalDomainElementBase<dim, codim,
      const PhysicalDomain<dim,codim>>
{
    using PhysicalDomainElementBase<dim, codim,
          const PhysicalDomain<dim,codim>>::PhysicalDomainElementBase;
};

template <int dim, int codim>
class PhysicalDomainElement
    : public PhysicalDomainElementBase<dim, codim,
      PhysicalDomain<dim,codim>>
{
    using PhysicalDomainElementBase<dim, codim,
          PhysicalDomain<dim,codim>>::PhysicalDomainElementBase;
};

IGA_NAMESPACE_CLOSE

#endif // PHYSICAL_DOMAIN_ELEMENT_H_


