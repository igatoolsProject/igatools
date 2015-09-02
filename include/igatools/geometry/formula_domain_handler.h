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

#include <igatools/geometry/domain_handler.h>
#include <igatools/geometry/formula_domain.h>

IGA_NAMESPACE_OPEN

/**
 *
 */
template<int dim, int codim>
class FormulaDomainHandler :
  public  DomainHandler<dim, codim>
{
private:
  using parent_t = DomainHandler<dim, codim>;
  using self_t = FormulaDomainHandler<dim, codim>;
protected:
  using typename parent_t::GridType;
public:
// using typename parent_t::DomainType;
  using DomainType =  const FormulaDomain<dim,codim>;
  using typename parent_t::ConstElementAccessor;

  using typename parent_t::topology_variant;
  using typename parent_t::eval_pts_variant;


//  using typename parent_t::Point;
//  using typename parent_t::Value;
//  using typename parent_t::Gradient;
//  using typename parent_t::ElementIterator;
//  using typename parent_t::ElementAccessor;
//  using parent_t::space_dim;

//  template <int order>
//  using Derivative = typename DomainType::template Derivative<order>;

  FormulaDomainHandler(std::shared_ptr<DomainType> domain);


  virtual ~FormulaDomainHandler() = default;


  void fill_cache(const topology_variant &sdim,
                  ConstElementAccessor &elem,
                  const int s_id) const override;


  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const DomainType &domain,
                        const self_t &domain_handler,
                        ConstElementAccessor &elem,
                        const int s_id)
      :
      domain_(domain),
      domain_handler_(domain_handler),
      elem_(elem),
      s_id_(s_id)
    {}


    template<int sdim>
    void operator()(const Topology<sdim> &sub_elem)
    {
      using _Point = typename ConstElementAccessor::_Point;


      auto &local_cache = domain_handler_.get_element_cache(elem_);
      auto &cache = local_cache->template get_sub_elem_cache<sdim>(s_id_);

      if (!cache.fill_none())
      {
        const auto &grid_pts = elem_.get_grid_element().template get_points<sdim>(s_id_);
        if (cache.template status_fill<domain_element::_Point>())
        {
          domain_.evaluate_0(grid_pts, cache.template get_data<_Point>());
          cache.template set_status_filled<domain_element::_Point>(true);
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
    const self_t     &domain_handler_;
    ConstElementAccessor &elem_;
    const int s_id_;
  };

  friend struct FillCacheDispatcher;
};

IGA_NAMESPACE_CLOSE

#endif
