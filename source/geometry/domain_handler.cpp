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

#include <igatools/geometry/domain_handler.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/functions/function.h>
#include <igatools/functions/function_handler.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
DomainHandler<dim_, codim_>::
DomainHandler(std::shared_ptr<DomainType> domain)
  :
  domain_(domain)
{
  Assert(domain_ != nullptr, ExcNullPtr());
}



template<int dim_, int codim_>
DomainHandler<dim_, codim_>::
~DomainHandler()
{}



template<int dim_, int codim_>
auto
DomainHandler<dim_, codim_>::
set_flags(const topology_variant &sdim,
          const Flags &flag) -> void
{
  using GridFlags = typename GridType::ElementHandler::Flags;
  GridFlags grid_flag = GridFlags::none;

  //point => grid::point
  grid_flag |= GridFlags::point;
  grid_handler_->set_flags(sdim, grid_flag);

  auto disp = SetFlagsDispatcher(flag, flags_);
  boost::apply_visitor(disp, sdim);

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
DomainHandler<dim_, codim_>::
init_cache(ConstElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  grid_handler_->init_cache(*(elem.grid_elem_), quad);


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
DomainHandler<dim_, codim_>::
fill_cache(const topology_variant &sdim,
           ConstElementAccessor &elem,
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

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/domain_handler.inst>

