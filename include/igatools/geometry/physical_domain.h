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

#ifndef __PHYSICAL_DOMAIN_H_
#define __PHYSICAL_DOMAIN_H_

#include <igatools/base/config.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/grid_cache_handler.h>

IGA_NAMESPACE_OPEN

namespace physical_domain_element
{
enum class Flags
{
    /** Fill nothing */
    none           =    0,

    /** Quadrature points on the element */
    points          =    1L << 1,

    /** Quadrature weigths on the element */
    w_measure       =    1L << 2
};
}

template <int, int, int, int> class Function;

template <int, int, class> class PhysicalDomainElementBase;
template <int, int> class PhysicalDomainElement;

/**
 * @brief The mapping is a deformation \f$ F : \hat\Omega \to \Omega\f$
 * which maps the reference domain \f$\hat\Omega \in \mathbb{R}^{dim}\f$ to the
 * physical domain \f$\Omega \in \mathbb{R}^{dim+codim}\f$.
 *
 * PhysicalDomain is the physical domain, wether of a function or a space.
 *
 * It is a function with special properties: it codim is 0 and the map is always the
 * identity.
 *
 * @todo we should thing about renaming mapping to physical domain
 *
 * @ingroup containers
 * @ingroup serializable
 *
 * @author pauletti 2014, 2015
 * @author M. Martinelli, 2015
 */
template<int dim_, int codim_ = 0>
class PhysicalDomain :
    public std::enable_shared_from_this<PhysicalDomain<dim_,codim_> >
{
private:
    using self_t = PhysicalDomain<dim_, codim_>;

public:
    static const int space_dim = dim_ + codim_;
    static const int dim = dim_;

    using GridType = const CartesianGrid<dim_>;
    using GridHandle = typename CartesianGrid<dim_>::ElementHandler;
    using FuncType =  Function<dim_, 0, dim_ + codim_, 1>;

    using ElementAccessor = PhysicalDomainElement<dim_, codim_>;
    using ElementIterator = GridIterator<ElementAccessor>;

    using List = typename GridType::List;
    using ListIt = typename GridType::ListIt;
    using Flags = physical_domain_element::Flags;

#if 0
public:
    /** Type for the given order derivatives of the
     *  the mapping. */
    template<int order>
    using Derivative = typename FuncType::template Derivative<order>;

    /** Type for the diferent order derivatives of the inverse of
     * the mapping
     */
    template<int order>
    using InvDerivative = Derivatives<space_dim, dim_, 1, order>;

    /** Type of the mapping evaluation point. */
    using Point = typename FuncType::Point;

    /** Type of the mapping return value. */
    using Value = typename FuncType::Value;

    /** Type of the mapping gradient. */
    using Gradient = typename FuncType::Gradient;

    /** Typedef for the mapping hessian. */
    using Hessian = typename FuncType::Hessian;

    using topology_variant = typename FuncType::topology_variant;
    using eval_pts_variant = typename FuncType::eval_pts_variant;
#endif
    using topology_variant = TopologyVariants<dim_>;
    using eval_pts_variant = SubElemPtrVariants<Quadrature,dim_>;
public:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    PhysicalDomain() = default;

    PhysicalDomain(std::shared_ptr<const GridType> grid);//,std::shared_ptr<const FuncType> F);

    ~PhysicalDomain();

//    static std::shared_ptr<self_t>  create(std::shared_ptr<const FuncType> F);
//
//    std::shared_ptr<self_t> get_handler() const
//    {
//        return std::make_shared<self_t>(self_t(this->F_->clone()));
//    }


    std::shared_ptr<const GridType> get_grid() const;

   // std::shared_ptr<const FuncType> get_function() const;

public:

    void set_flags(const topology_variant &sdim,
                   const Flags &flag);

    void init_cache(ElementAccessor &elem,
                    const eval_pts_variant &quad) const;

    void init_cache(ElementIterator &elem,
                    const eval_pts_variant &quad) const
    {
        this->init_cache(*elem, quad);
    }

    void fill_cache(const topology_variant &sdim,
                    ElementAccessor &elem,
                    const int s_id) const;

    void fill_cache(const topology_variant &sdim,
                    ElementIterator &elem,
                    const int s_id) const
    {
        this->fill_cache(sdim, *elem, s_id);
    }

    ElementIterator begin(const PropId &property = ElementProperties::active);

    ElementIterator end(const PropId &property = ElementProperties::active);

#if 0
private:
    struct ResetDispatcher : boost::static_visitor<void>
    {
        ResetDispatcher(const ValueFlags flag_in,
                        SafeSTLArray<ValueFlags, dim_ + 1> &flags)
            :
            flag_in_(flag_in),
            flags_(flags)
        {}

        template<int sub_elem_dim>
        void operator()(const Quadrature<sub_elem_dim> &quad)
        {
            flags_[sub_elem_dim] = flag_in_;
        }

        const ValueFlags flag_in_;
        SafeSTLArray<ValueFlags, dim_ + 1> &flags_;
    };


    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        FillCacheDispatcher(const FuncType &F,
                            ElementAccessor &domain_elem,
                            const int sub_elem_id)
            :
            F_(F),
            domain_elem_(domain_elem),
            s_id_(sub_elem_id)

        {}

        template<int sdim>
        void operator()(const Topology<sdim> &sub_elem)
        {
            const int s_id=s_id_;
            auto &elem = domain_elem_;

            // TODO (pauletti, Nov 6, 2014): provide a lighter function for this
            const auto n_points = elem.get_grid_elem()->template get_quad<sdim>()->
            get_num_points();

            auto &cache = elem.local_cache_->template get_sub_elem_cache<sdim>(s_id);

            if (cache.template status_fill<_Measure>())
            {
                auto &s_elem = UnitElement<dim_>::template get_elem<sdim>(s_id);

                const auto &DF = elem.template get_values<_Gradient, sdim>(s_id);
                typename MapFunction<sdim, space_dim>::Gradient DF1;

                auto &measures = cache.template get_data<_Measure>();
                for (int pt = 0 ; pt < n_points; ++pt)
                {
                    for (int l=0; l<sdim; ++l)
                        DF1[l] = DF[pt][s_elem.active_directions[l]];

                    measures[pt] = fabs(determinant<sdim,space_dim>(DF1));
                }
                cache.template set_status_filled<_Measure>(true);
            }

            if (cache.template status_fill<_W_Measure>())
            {
                const auto &w = elem.grid_elem_->template get_weights<sdim>(s_id);

                const auto &measures = cache.template get_data<_Measure>();

                auto &w_measures = cache.template get_data<_W_Measure>();

                for (int pt = 0 ; pt < n_points; ++pt)
                    w_measures[pt] = w[pt] * measures[pt];

                cache.template set_status_filled<_W_Measure>(true);
            }

            if (cache.template status_fill<_InvGradient>())
            {
                // TODO (pauletti, Nov 23, 2014): if also fill measure this could be done here
                const auto &DF = elem.template get_values<_Gradient, sdim>(s_id);
                auto &D_invF = cache.template get_data<_InvGradient>();
                Real det;
                for (int pt = 0 ; pt < n_points; ++pt)
                    D_invF[pt] = inverse(DF[pt], det);

                cache.template set_status_filled<_InvGradient>(true);
            }

            if (cache.template status_fill<_InvHessian>())
            {
                //        const auto &D1_F = elem.template get_values<_Gradient, k>(j);
                const auto &D2_F = elem.template get_values<_Hessian, sdim>(s_id);
                const auto &D1_invF = cache.template get_data<_InvGradient>();
                auto &D2_invF       = cache.template get_data<_InvHessian>();

                for (int pt = 0 ; pt < n_points; ++pt)
                    for (int u=0; u<dim_; ++u)
                    {
                        const auto tmp_u = action(D2_F[pt], D1_invF[pt][u]);
                        for (int v=0; v<dim_; ++v)
                        {
                            const auto tmp_u_v = action(tmp_u, D1_invF[pt][v]);
                            D2_invF[pt][u][v] = - action(D1_invF[pt], tmp_u_v);
                        }
                    }

                cache.template set_status_filled<_InvHessian>(true);
            }

            if (cache.template status_fill<_BoundaryNormal>())
            {
                Assert(dim_ == sdim+1, ExcNotImplemented());
                const auto &D1_invF = cache.template get_data<_InvGradient>();
                const auto n_hat  = F_.get_grid()->template get_boundary_normals<sdim>(s_id)[0];
                auto &bndry_normal = cache.template get_data<_BoundaryNormal>();

                for (int pt = 0; pt < n_points; ++pt)
                {
                    const auto D1_invF_t = co_tensor(transpose(D1_invF[pt]));
                    bndry_normal[pt] = action(D1_invF_t, n_hat);
                    bndry_normal[pt] /= bndry_normal[pt].norm();
                }

                cache.template set_status_filled<_BoundaryNormal>(true);
            }

            if (cache.template status_fill<_OuterNormal>())
            {
                Assert(sdim == dim_, ExcNotImplemented());
                Assert(codim_ == 1, ExcNotImplemented());

                const auto &DF = elem.template get_values<_Gradient, sdim>(s_id);
                auto &outer_normal = cache.template get_data<_OuterNormal>();

                for (int pt = 0; pt < n_points; ++pt)
                {
                    outer_normal[pt] = cross_product<dim_, codim_>(DF[pt]);
                    outer_normal[pt] /= outer_normal[pt].norm();
                }

                cache.template set_status_filled<_OuterNormal>(true);
            }


            if (cache.template status_fill<_Curvature>())
            {
                Assert(sdim == dim_, ExcNotImplemented());
                Assert(codim_ == 1, ExcNotImplemented());

                const auto H = elem.compute_second_fundamental_form();
                const auto G_inv = elem.compute_inv_first_fundamental_form();

                auto &curvatures = cache.template get_data<_Curvature>();

                for (int pt = 0; pt < n_points; ++pt)
                {
                    //          const MetricTensor B = compose(H[pt], G_inv[pt]);
                    const auto B = compose(H[pt], G_inv[pt]);
                    const auto A = unroll_to_matrix(B);
                    curvatures[pt] = A.eigen_values();
                }

                cache.template set_status_filled<_Curvature>(true);
            }
            cache.set_filled(true);
        }

        const FuncType &F_;
        ElementAccessor &domain_elem_;
        int s_id_;
    };


    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        InitCacheDispatcher(const FuncType &F,
                            ElementAccessor &domain_elem,
                            const SafeSTLArray<ValueFlags, dim_ + 1> &flags)
            :
            F_(F),
            domain_elem_(domain_elem),
            flags_(flags)
        {}


        template<int k>
        void operator()(const Topology<k> &sub_elem)
        {
            auto &cache = domain_elem_.local_cache_;
            for (auto &s_id: UnitElement<dim_>::template elems_ids<k>())
            {
                auto &s_cache = cache->template get_sub_elem_cache<k>(s_id);
                const auto n_points = F_.template get_num_points<k>();
                s_cache.resize(flags_[k], n_points);
            }
        }

        const FuncType &F_;
        ElementAccessor &domain_elem_;
        const SafeSTLArray<ValueFlags, dim_ + 1> &flags_;
    };
#endif

private:
    std::shared_ptr<const GridType> grid_;
    std::shared_ptr<GridHandle> grid_handler_;
    //std::shared_ptr<const FuncType> F_;

    SafeSTLArray<Flags, dim_ + 1> flags_;

    friend ElementAccessor;

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
        ar.template register_type<IgFunction<dim_,0,dim_+codim_,1> >();
        ar &boost::serialization::make_nvp("F_",F_);
        ar &boost::serialization::make_nvp("flags_",flags_);
    }
    ///@}
#endif
};

IGA_NAMESPACE_CLOSE

#endif

