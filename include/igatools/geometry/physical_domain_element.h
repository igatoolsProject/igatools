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

#include <igatools/utils/safe_stl_array.h>
#include <igatools/geometry/physical_domain.h>
#include <igatools/functions/function_element.h>

IGA_NAMESPACE_OPEN

/**
 *
 * @ingroup elements
 */
template<int dim_, int codim_ = 0>
class PhysicalDomainElement
{
private:
    using self_t  = PhysicalDomainElement<dim_, codim_>;
    using parent_t = FunctionElement<dim_, 0, dim_+codim_>;
    using PhysDom = PhysicalDomain<dim_, codim_>;
    using Func = MapFunction<dim_, codim_>;   

public:
    using ContainerType = PhysDom;
    static const int dim = dim_;
    static const int codim = codim_;
    static const int space_dim = dim_+codim_;


    /** @name Constructors */
    ///@{
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    PhysicalDomainElement() = default;

    /**
     * Construct an accessor pointing to the element with
     * flat index @p elem_index of the Function @p func.
     */
    PhysicalDomainElement(const std::shared_ptr<const PhysDom> func,
                          const Index elem_index);

    /**
     * Copy constructor.
     * It can be used with different copy policies
     * (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a
     * classic copy constructor)
     * uses the deep copy.
     */
    PhysicalDomainElement(const self_t &elem,
                          const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    PhysicalDomainElement(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~PhysicalDomainElement() = default;
    ///@}

    template<int order>
    using InvDerivative = typename Map::template InvDerivative<order>;

    template <int order>
    using Derivative = typename Map::template Derivative<order>;

private:
    template <class ValueType, int topology_dim = dim>
    auto &get_values_from_cache(const int topology_id = 0) const
    {
        Assert(local_cache_ != nullptr,ExcNullPtr());
        const auto &cache = local_cache_->template
        		get_sub_elem_cache<topology_dim>(topology_id);
        return cache.template get_data<ValueType>();
    }

public:
    template<int k>
    ValueVector<Real> const &get_measures(const int j) const
    {
        return get_values_from_cache<_Measure,k>(j);
    }

    template<int k>
    ValueVector<Real> const &get_w_measures(const int j) const
    {
        return get_values_from_cache<_W_Measure,k>(j);
    }

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


private:

    using CType = boost::fusion::map<
                  boost::fusion::pair<       _Measure,DataWithFlagStatus<ValueVector<Real>>>,
                  boost::fusion::pair<     _W_Measure,DataWithFlagStatus<ValueVector<Real>>>,
                  boost::fusion::pair<   _InvGradient,DataWithFlagStatus<ValueVector<InvDerivative<1>>>>,
                  boost::fusion::pair<    _InvHessian,DataWithFlagStatus<ValueVector<InvDerivative<2>>>>,
                  boost::fusion::pair<_BoundaryNormal,DataWithFlagStatus<ValueVector<Points<space_dim>>>>,
                  boost::fusion::pair<   _OuterNormal,DataWithFlagStatus<ValueVector<Points<space_dim>>>>,
                  boost::fusion::pair<     _Curvature,DataWithFlagStatus<ValueVector<SafeSTLVector<Real>>>>
                  >;


    /**
     * Returns the flags that are valid to be used with this class.
     *
     * @note The valid flags are defined to be the ones that can be inferred from the ValueType(s)
     * used as key of the boost::fusion::map in CType.
     */
    static ValueFlags get_valid_flags()
    {
        return cacheutils::get_valid_flags_from_cache_type(CType());
    }

    using Cache = FuncValuesCache<dim,CType>;


public:
    using CacheType = AllSubElementsCache<Cache>;

private:

    using FuncElem = FunctionElement<dim_, 0, dim_+codim_>;

    std::shared_ptr<FuncElem> func_elem_;

    std::shared_ptr<CacheType> local_cache_;


    template <class Accessor> friend class CartesianGridIteratorBase;
    friend class PhysicalDomain<dim, codim>;

    /**
     * Creates a new object performing a deep copy of the current object using the PhysicalDomainElement
     * copy constructor.
     */
    std::shared_ptr<PhysicalDomainElement<dim_,codim_> > clone() const;

};


IGA_NAMESPACE_CLOSE

#endif // PHYSICAL_DOMAIN_ELEMENT_H_



