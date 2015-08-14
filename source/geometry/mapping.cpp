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

#include <igatools/base/tensor.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN


namespace
{

ValueFlags
mapping_to_function_flags(const ValueFlags &flags)
{
	/*
    ValueFlags valid_func_flags = ValueFlags::gradient |
                                  ValueFlags::hessian |
                                  ValueFlags::divergence |
                                  ValueFlags::point;
    ValueFlags transfer_flags = ValueFlags::measure |
                                ValueFlags::w_measure |
                                ValueFlags::boundary_normal |
                                valid_func_flags;
//*/

	ValueFlags transfer_flags = ValueFlags::gradient |
			                    ValueFlags::hessian;
    ValueFlags f_flags = flags & transfer_flags;

    if (contains(flags, ValueFlags::point))
        f_flags |= ValueFlags::value;


    if (contains(flags, ValueFlags::measure) ||
        contains(flags, ValueFlags::w_measure) ||
        contains(flags, ValueFlags::inv_gradient) ||
        contains(flags, ValueFlags::outer_normal))
        f_flags |=  ValueFlags::gradient;

    if (contains(flags, ValueFlags::inv_hessian) ||
        contains(flags, ValueFlags::curvature))
        f_flags |=  ValueFlags::gradient | ValueFlags::hessian;

    return f_flags;
}
};



template<int dim_, int codim_>
Mapping<dim_, codim_>::
Mapping(std::shared_ptr<const FuncType> F)
    :
    F_(F->clone())
{}



template<int dim_, int codim_>
Mapping<dim_, codim_>::
~Mapping()
{}



template<int dim_, int codim_>
auto
Mapping<dim_, codim_>::
create(std::shared_ptr<const FuncType> F)-> std::shared_ptr<self_t>
{
    return std::shared_ptr<self_t>(new self_t(F));
}


template<int dim_, int codim_>
template <int k>
auto
Mapping<dim_, codim_>::
reset(const ValueFlags flags, const Quadrature<k> &eval_pts) -> void
{
    const auto valid_flags = ElementAccessor::get_valid_flags();
    auto m_flags = flags & valid_flags;

    if (contains(flags, ValueFlags::boundary_normal) ||
    contains(flags, ValueFlags::curvature))
        m_flags |= ValueFlags::inv_gradient;


    if (contains(flags, ValueFlags::curvature))
        m_flags |= ValueFlags::outer_normal;


    if (contains(flags, ValueFlags::w_measure))
        m_flags |= ValueFlags::measure;


    std::const_pointer_cast<FuncType>(F_)->reset(mapping_to_function_flags(m_flags), eval_pts);
    flags_[k] = m_flags;
}



template<int dim_, int codim_>
template <int k>
auto
Mapping<dim_, codim_>::
fill_cache(ElementAccessor &elem, const int j) const -> void
{
    F_->template fill_cache(elem.get_func_element(), Topology<k>(),j);

    // TODO (pauletti, Nov 6, 2014): provide a lighter function for this
    const auto n_points = F_->template get_num_points<k>();

    auto &cache = elem.local_cache_->template get_sub_elem_cache<k>(j);
    const auto &func_elem = elem.get_func_element();

    if (cache.template status_fill<_Point>())
    {
        cache.template get_data<_Point>() =
        		func_elem.template get_values<_Value,k>(j);
    }

    if (cache.template status_fill<_Gradient>())
    {
        cache.template get_data<_Gradient>() =
        		func_elem.template get_values<_Gradient,k>(j);
    }

    if (cache.template status_fill<_Hessian>())
    {
        cache.template get_data<_Hessian>() =
        		func_elem.template get_values<_Hessian,k>(j);
    }

    if (cache.template status_fill<_Measure>())
    {
        auto &k_elem = UnitElement<dim_>::template get_elem<k>(j);

        const auto &DF = func_elem.template get_values<_Gradient, k>(j);
//        typename MapFunction<k, space_dim>::Gradient DF1;
        typename MapFunction_new<k, space_dim - k>::Gradient DF1;

        auto &measures = cache.template get_data<_Measure>();
        for (int pt = 0 ; pt < n_points; ++pt)
        {
            for (int l=0; l<k; ++l)
                DF1[l] = DF[pt][k_elem.active_directions[l]];

            measures[pt] = fabs(determinant<k,space_dim>(DF1));
        }
        cache.template set_status_filled<_Measure>(true);
    }

    if (cache.template status_fill<_W_Measure>())
    {
        const auto &w = func_elem.get_grid_element().template get_w_measures<k>(j);

        const auto &measures = cache.template get_data<_Measure>();

        auto &w_measures = cache.template get_data<_W_Measure>();

        for (int pt = 0 ; pt < n_points; ++pt)
            w_measures[pt] = w[pt] * measures[pt];

        cache.template set_status_filled<_W_Measure>(true);
    }

    if (cache.template status_fill<_InvGradient>())
    {
        // TODO (pauletti, Nov 23, 2014): if also fill measure this could be done here
        const auto &DF = func_elem.template get_values<_Gradient, k>(j);
        auto &D_invF = cache.template get_data<_InvGradient>();
        Real det;
        for (int pt = 0 ; pt < n_points; ++pt)
            D_invF[pt] = inverse(DF[pt], det);

        cache.template set_status_filled<_InvGradient>(true);
    }

    if (cache.template status_fill<_InvHessian>())
    {
//        const auto &D1_F = elem.template get_values<_Gradient, k>(j);
        const auto &D2_F = func_elem.template get_values<_Hessian, k>(j);
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
        Assert(dim_ == k+1, ExcNotImplemented());
        const auto &D1_invF = cache.template get_data<_InvGradient>();
        const auto n_hat  = this->get_grid()->template get_boundary_normals<k>(j)[0];
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
        Assert(k == dim_, ExcNotImplemented());
        Assert(codim_ == 1, ExcNotImplemented());

        const auto &DF = func_elem.template get_values<_Gradient, k>(j);
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
        Assert(k == dim_, ExcNotImplemented());
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



template<int dim_, int codim_>
template <int k>
auto
Mapping<dim_, codim_>::
init_cache(ElementAccessor &elem) const -> void
{
    F_->init_cache(elem.get_func_element(), Topology<k>());

    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::CacheType;
        cache = shared_ptr<Cache>(new Cache);
    }

    for (auto &s_id: UnitElement<dim_>::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_sub_elem_cache<k>(s_id);
        const auto n_points = F_->template get_num_points<k>();
        s_cache.resize(flags_[k], n_points);
    }

}


//    if (flag_.fill_inv_hessians())
//    {
//        const auto &D1_F = elem.get_gradients();
//        const auto &D2_F = elem.get_hessians();
//        const auto &D1_invF = std::get<1>(cache->inv_derivatives_);
//        auto &D2_invF = std::get<2>(cache->inv_derivatives_);
//
//        for (int i=0; i<n_points; ++i)
//            for (int u=0; u<dim_; ++u)
//            {
//                const auto tmp_u = action(D2_F[i], D1_invF[i][u]);
//                for (int v=0; v<dim_; ++v)
//                {
//                    const auto tmp_u_v = action(tmp_u, D1_invF[i][v]);
//                    D2_invF[i][u][v] = - action(D1_invF[i], tmp_u_v);
//                }
//            }
//    }
//}


template<int dim_, int codim_>
auto
Mapping<dim_, codim_>::
get_grid() const -> std::shared_ptr<const CartesianGrid<dim_> >
{
    return F_->get_grid();
}

template<int dim_, int codim_>
auto
Mapping<dim_, codim_>::
get_function() const -> std::shared_ptr<const FuncType>
{
    return F_;
}

template<int dim_, int codim_>
auto
Mapping<dim_, codim_>::
create_element(const Index flat_index) const -> std::shared_ptr<ElementAccessor>
{
    auto elem = std::make_shared<ElementAccessor>(this->get_function(),flat_index);
    Assert(elem != nullptr, ExcNullPtr());

    return elem;
}


template<int dim_, int codim_>
auto
Mapping<dim_, codim_>::
begin() const -> ElementIterator
{
    return ElementIterator(this->create_element(0),ElementProperties::none);
}

template<int dim_, int codim_>
auto
Mapping<dim_, codim_>::
end() -> ElementIterator
{
    return ElementIterator(this->create_element(IteratorState::pass_the_end),ElementProperties::none);
}


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/mapping.inst>

