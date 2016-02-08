//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#ifndef __FORMULA_FUNCTION_H_
#define __FORMULA_FUNCTION_H_

#include <igatools/functions/function.h>
#include <igatools/functions/function_element.h>

IGA_NAMESPACE_OPEN

template <int, int, int, int> class FormulaFunctionHandler;

/**
 *
 */
template<int dim, int codim=0, int range = 1, int rank = 1>
class FormulaFunction :
  public Function<dim, codim, range, rank>
{
private:
  using parent_t = Function<dim, codim, range, rank>;
  using self_t = FormulaFunction<dim, codim, range, rank>;
protected:
  using typename parent_t::DomainType;
  using Handler = FormulaFunctionHandler<dim, codim, range, rank>;
public:
  using typename parent_t::Value;
  using typename parent_t::Point;

  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  FormulaFunction(const SharedPtrConstnessHandler<DomainType> &domain);

// FormulaFunction(const self_t &func);

  virtual ~FormulaFunction() = default;

  std::unique_ptr<typename parent_t::Handler>
  create_cache_handler() const override final;

public:

  virtual void evaluate_0(const ValueVector<Point> &points,
                          ValueVector<Value> &values) const = 0;

  virtual void evaluate_1(const ValueVector<Point> &points,
                          ValueVector<Derivative<1>> &values) const = 0;

  virtual void evaluate_2(const ValueVector<Point> &points,
                          ValueVector<Derivative<2>> &values) const = 0;



  virtual void print_info(LogStream &out) const override = 0;

#if 0
  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const int sub_elem_id,const self_t &function,ElementAccessor &elem)
      :
      sub_elem_id_(sub_elem_id),
      function_(function),
      elem_(elem)
    {}


    template<int sub_elem_dim>
    void operator()(const Topology<sub_elem_dim> &sub_elem)
    {
      auto &local_cache = elem_.get_cache();
      auto &cache = local_cache->template get_sub_elem_cache<sub_elem_dim>(sub_elem_id_);

      if (!cache.fill_none())
      {
        auto &cache_pts = cache.template get_data<_Point>();
        cache_pts = function_.map_elem_->template get_values<_Value, sub_elem_dim>(sub_elem_id_);
        if (cache.template status_fill<_Value>())
        {
          function_.evaluate_0(cache_pts, cache.template get_data<_Value>());
          cache.template set_status_filled<_Value>(true);
        }
        if (cache.template status_fill<_Gradient>())
        {
          function_.evaluate_1(cache_pts, cache.template get_data<_Gradient>());
          cache.template set_status_filled<_Gradient>(true);
        }
        if (cache.template status_fill<_Hessian>())
        {
          function_.evaluate_2(cache_pts, cache.template get_data<_Hessian>());
          cache.template set_status_filled<_Hessian>(true);
        }
        if (cache.template status_fill<_Divergence>())
          Assert(false,ExcNotImplemented());
      }

      cache.set_filled(true);
    }

    const int sub_elem_id_;
    const self_t &function_;
    ElementAccessor &elem_;
  };

  friend struct FillCacheDispatcher;
#endif

};

IGA_NAMESPACE_CLOSE

#endif

