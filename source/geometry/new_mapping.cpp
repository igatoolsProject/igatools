//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/mapping_element.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim, int codim>
auto
NewMapping<dim, codim>::
get_cache(ElementAccessor &elem)
-> std::shared_ptr<typename ElementAccessor::CacheType> &
{
    return elem.local_cache_;
}



template<int dim, int codim>
NewMapping<dim, codim>::
NewMapping(std::shared_ptr<FuncType> F)
    :
    F_(F)
{}



template<int dim, int codim>
NewMapping<dim, codim>::
~NewMapping()
{}



template<int dim, int codim>
auto
NewMapping<dim, codim>::
create(std::shared_ptr<FuncType> F)-> std::shared_ptr<self_t>
{
    return std::shared_ptr<self_t>(new self_t(F));
}


template<int dim, int codim>
template <int k>
auto
NewMapping<dim, codim>::
reset(const NewValueFlags flag, const Quadrature<k> &quad) -> void
{
    F_->reset(MappingFlags::to_function_flags(flag), quad);
    flags_[k] = flag;
}



template<int dim, int codim>
template <int k>
auto
NewMapping<dim, codim>::
fill_cache(ElementAccessor &elem, const int j) -> void
{
    F_->template fill_cache(elem, j, Int<k>());

    // TODO (pauletti, Nov 6, 2014): provide a lighter function for this
    const auto n_points = F_->template get_num_points<k>();

    auto &cache = elem.local_cache_->template get_value_cache<k>(j);
    auto &flags = cache.flags_handler_;

    if (flags.fill_measures())
    {
        auto &k_elem = UnitElement<dim>::template get_elem<k>(j);

        const auto &DF = elem.template get_values<1, k>(j);
        typename MapFunction<k, space_dim>::Gradient DF1;


        for (int i=0; i<n_points; ++i)
        {
            for (int l=0; l<k; ++l)
                DF1[l] = DF[i][k_elem.active_directions[l]];

            cache.measures_[i] = fabs(determinant<k,space_dim>(DF1));
        }
    }

    if (flags.fill_w_measures())
    {
        const auto &meas = cache.measures_;
        const auto &w = elem.CartesianGridElement<dim>::template get_w_measures<k>(j);
        for (int i=0; i<n_points; ++i)
            cache.w_measures_[i] = w[i] * meas[i];
    }

    if (flags.fill_inv_gradients())
    {
    	// TODO (pauletti, Nov 23, 2014): if also fill measure this could be done here
        const auto &DF = elem.template get_values<1, k>(j);
        auto &D_invF = cache.template get_inv_values<1>();
        Real det;
        for (int i=0; i<n_points; ++i)
        	D_invF[i] = inverse(DF[i], det);
    }

    if (flags.fill_inv_hessians())
    {
        const auto &D1_F = elem.template get_values<1, k>(j);
        const auto &D2_F = elem.template get_values<2, k>(j);
        const auto &D1_invF = cache.template get_inv_values<1>();
        auto &D2_invF       = cache.template get_inv_values<2>();

        for (int i=0; i<n_points; ++i)
            for (int u=0; u<dim; ++u)
            {
                const auto tmp_u = action(D2_F[i], D1_invF[i][u]);
                for (int v=0; v<dim; ++v)
                {
                    const auto tmp_u_v = action(tmp_u, D1_invF[i][v]);
                    D2_invF[i][u][v] = - action(D1_invF[i], tmp_u_v);
                }
            }
    }

    cache.set_filled(true);
}



template<int dim, int codim>
template <int k>
auto
NewMapping<dim, codim>::
init_cache(ElementAccessor &elem) -> void
{
    F_->init_cache(elem, Int<k>());

    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    for (auto &s_id: UnitElement<dim>::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_value_cache<k>(s_id);
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
//            for (int u=0; u<dim; ++u)
//            {
//                const auto tmp_u = action(D2_F[i], D1_invF[i][u]);
//                for (int v=0; v<dim; ++v)
//                {
//                    const auto tmp_u_v = action(tmp_u, D1_invF[i][v]);
//                    D2_invF[i][u][v] = - action(D1_invF[i], tmp_u_v);
//                }
//            }
//    }
//}


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/new_mapping.inst>

