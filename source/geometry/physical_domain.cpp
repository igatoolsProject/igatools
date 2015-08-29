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

#include <igatools/geometry/physical_domain.h>
#include <igatools/geometry/physical_domain_element.h>
#include <igatools/functions/function.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN


//namespace
//{
//
//ValueFlags
//mapping_to_function_flags(const ValueFlags &flags)
//{
//    ValueFlags valid_func_flags = ValueFlags::value |
//                                  ValueFlags::gradient |
//                                  ValueFlags::hessian |
//                                  ValueFlags::divergence |
//                                  ValueFlags::point;
//
//    ValueFlags transfer_flags = ValueFlags::measure |
//                                ValueFlags::w_measure |
//                                ValueFlags::boundary_normal |
//                                valid_func_flags;
//
//
//    ValueFlags f_flags = flags & transfer_flags;
//
//    if (contains(flags, ValueFlags::measure) ||
//        contains(flags, ValueFlags::w_measure) ||
//        contains(flags, ValueFlags::inv_gradient) ||
//        contains(flags, ValueFlags::outer_normal))
//        f_flags |=  ValueFlags::gradient;
//
//    if (contains(flags, ValueFlags::inv_hessian) ||
//        contains(flags, ValueFlags::curvature))
//        f_flags |=  ValueFlags::gradient | ValueFlags::hessian;
//
//    return f_flags;
//}
//};



template<int dim_, int codim_>
PhysicalDomain<dim_, codim_>::
PhysicalDomain(std::shared_ptr<const GridType> grid,
               std::shared_ptr<const FuncType> F)
    :
    grid_(grid),
    grid_handler_(grid->create_cache_handler()),
    func_(F)
{}



template<int dim_, int codim_>
PhysicalDomain<dim_, codim_>::
~PhysicalDomain()
{}



//template<int dim_, int codim_>
//auto
//PhysicalDomain<dim_, codim_>::
//create(std::shared_ptr<FuncType> F)-> std::shared_ptr<self_t>
//{
//    return std::shared_ptr<self_t>(new self_t(F));
//}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
set_flags(const topology_variant &sdim,
          const Flags &flag) -> void
{
#if 0
    const auto valid_flags = ElementAccessor::get_valid_flags();
    auto m_flags = flags & valid_flags;

    if (contains(flags, ValueFlags::boundary_normal) ||
    contains(flags, ValueFlags::curvature))
        m_flags |= ValueFlags::inv_gradient;

    if (contains(flags, ValueFlags::curvature))
        m_flags |= ValueFlags::outer_normal;

    if (contains(flags, ValueFlags::w_measure))
        m_flags |= ValueFlags::measure;

    F_->reset(mapping_to_function_flags(m_flags), eval_pts);

    auto reset_dispatcher = ResetDispatcher(m_flags, flags_);
    boost::apply_visitor(reset_dispatcher, eval_pts);
#endif
}




template<int dim_, int codim_>
void
PhysicalDomain<dim_, codim_>::
init_cache(ElementAccessor &elem,
           const eval_pts_variant &quad) const
{
#if 0
    F_->init_cache(elem, k);

    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::CacheType;
        cache = shared_ptr<Cache>(new Cache);
    }

    auto init_cache_dispatcher = InitCacheDispatcher(*F_, elem, flags_);
    boost::apply_visitor(init_cache_dispatcher, k);
#endif
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
fill_cache(const topology_variant &sdim,
           ElementAccessor &elem,
           const int s_id) const-> void
{
#if 0
    F_->template fill_cache(elem, k, j);
    auto fill_cache_dispatcher =FillCacheDispatcher(*F_, elem, j);
    boost::apply_visitor(fill_cache_dispatcher, k);
#endif
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
PhysicalDomain<dim_, codim_>::
get_grid() const -> std::shared_ptr<const CartesianGrid<dim_> >
{
    return grid_;
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
get_function() const -> std::shared_ptr<const FuncType>
{
    return func_;
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
create_element(const ListIt &index, const PropId &prop) const
-> std::shared_ptr<ConstElementAccessor>
{
    using Elem = ConstElementAccessor;
    auto elem = std::make_shared<Elem>(this->shared_from_this(), index, prop);
    Assert(elem != nullptr,ExcNullPtr());

    return elem;
}



template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
begin(const PropId &prop) -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
    grid_->get_element_property(prop).begin(),
    prop);
}

template<int dim_, int codim_>
auto
PhysicalDomain<dim_, codim_>::
end(const PropId &prop) -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
    grid_->get_element_property(prop).end(),
    prop);
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/physical_domain.inst>

