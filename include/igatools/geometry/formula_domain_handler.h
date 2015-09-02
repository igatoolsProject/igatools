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

#ifndef __FORMULA_DOMAIN_HANDLER_H_
#define __FORMULA_DOMAIN_HANDLER_H_

#include <igatools/geometry/physical_domain_cache_handler.h>
#include <igatools/geometry/formula_domain.h>

IGA_NAMESPACE_OPEN

/**
 *
 */
template<int dim, int codim>
class FormulaDomainHandler :
  public  PhysicalDomainElementHandler<dim, codim>
{
private:
  using parent_t = PhysicalDomainElementHandler<dim, codim>;
  using self_t = FormulaDomainHandler<dim, codim>;
protected:
  using typename parent_t::GridType;
public:
  using typename parent_t::DomainType;

  using typename parent_t::topology_variant;
  using typename parent_t::eval_pts_variant;
//  using typename parent_t::Point;
//  using typename parent_t::Value;
//  using typename parent_t::Gradient;
//  using typename parent_t::ElementIterator;
//  using typename parent_t::ElementAccessor;
//  using parent_t::space_dim;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  FormulaDomainHandler(std::shared_ptr<DomainType> domain);

  //FormulaDomainHandler(const self_t &func);

  virtual ~FormulaDomainHandler() = default;

//  void set_flags(const topology_variant &sdim,
//                          const Flags &flag) override;

//  void init_cache(ElementAccessor &elem,
                           const eval_pts_variant &quad) const override;

  void fill_cache(const topology_variant &sdim,
                           ElementAccessor &elem,
                           const int s_id) const override;


  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const DomainType &domain,
    ConstElementAccessor &elem,
    const int s_id)
      :
    	  domain_(domain),
		  elem_(elem),
		  s_id_(sub_elem_id)
    {}


    template<int sdim>
    void operator()(const Topology<sdim> &sub_elem)
    {
    	using _Point = typename ConstElementAccessor::_Point;

      auto &local_cache = elem_.get_cache();
      auto &cache = local_cache->template get_sub_elem_cache<sdim>(s_id_);

      if (!cache.fill_none())
      {
    	  auto &grid_pts = elem_.get_grid_element()->get_points(s_id_);
    	  if (cache.template status_fill<_Value>())
    	  {
    		  domain_.evaluate_0(grid_pts, cache.template get_data<_Value>());
    		  cache.template set_status_filled<_Value>(true);
    	  }
//        if (cache.template status_fill<_Gradient>())
//        {
//          function_.evaluate_1(cache_pts, cache.template get_data<_Gradient>());
//          cache.template set_status_filled<_Gradient>(true);
//        }
//        if (cache.template status_fill<_Hessian>())
//        {
//          function_.evaluate_2(cache_pts, cache.template get_data<_Hessian>());
//          cache.template set_status_filled<_Hessian>(true);
//        }
//        if (cache.template status_fill<_Divergence>())
//          Assert(false,ExcNotImplemented());
      }

      cache.set_filled(true);
    }

    const DomainType &domain_;
    ConstElementAccessor &elem_;
    const int s_id_;
  };

  friend struct FillCacheDispatcher;
#endif

};

IGA_NAMESPACE_CLOSE

#endif
#endif
