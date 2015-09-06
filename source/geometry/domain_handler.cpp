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
//#include <igatools/functions/function.h>
//#include <igatools/functions/function_handler.h>

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
DomainHandler<dim_, codim_>::
DomainHandler(std::shared_ptr<DomainType> domain)
  :
  domain_(domain),
  grid_handler_(domain->get_grid()->create_cache_handler())
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
  GridFlags  grid_flag = GridFlags::none;
  CacheFlags dom_flag = CacheFlags::none;

  SafeSTLVector<Flags> all_flags ={Flags::point, Flags::measure, Flags::w_measure};
  for (auto &fl : all_flags)
    if (contains(flag, fl))
    {
      grid_flag |= domain_element::activate::grid[fl];
      dom_flag  |= domain_element::activate::domain[fl];
    }

  grid_handler_->set_flags(sdim, grid_flag);

  auto disp = SetFlagsDispatcher(dom_flag, flags_);
  boost::apply_visitor(disp, sdim);
}



template<int dim_, int codim_>
void
DomainHandler<dim_, codim_>::
init_cache(ConstElementAccessor &elem,
           const eval_pts_variant &quad) const
{
  grid_handler_->init_cache(*(elem.grid_elem_), quad);

  auto &cache = elem.local_cache_;
  if (cache == nullptr)
  {
    using Cache = typename ElementAccessor::CacheType;
    cache = std::make_shared<Cache>();
  }

  auto disp = InitCacheDispatcher(this, elem, flags_);
  boost::apply_visitor(disp, quad);

}



template<int dim_, int codim_>
auto
DomainHandler<dim_, codim_>::
fill_cache(const topology_variant &sdim,
           ConstElementAccessor &elem,
           const int s_id) const-> void
{
  grid_handler_->fill_cache(sdim, *(elem.grid_elem_), s_id);

  auto disp = FillCacheDispatcher(elem, s_id);
  boost::apply_visitor(disp, sdim);

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

