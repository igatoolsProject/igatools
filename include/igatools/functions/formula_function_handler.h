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

#ifndef __FORMULA_FUNCTION_HANDLER_H_
#define __FORMULA_FUNCTION_HANDLER_H_

#include <igatools/functions/function_handler.h>
#include <igatools/functions/formula_function.h>

IGA_NAMESPACE_OPEN

/**
 * @ingroup handlers
 */
template<int dim, int codim=0, int range = 1, int rank = 1>
class FormulaFunctionHandler :
  public FunctionHandler<dim, codim, range, rank>
{
private:
  using parent_t = FunctionHandler<dim, codim, range, rank>;
  using self_t = FormulaFunctionHandler<dim, codim, range, rank>;
protected:
  using typename parent_t::DomainType;
  ///using typename Handler = FormulaFunctionHandler<dim, codim, range, rank>;
public:
  using FuncType = const FormulaFunction<dim, codim, range, rank>;
  using typename parent_t::ElementAccessor;
  using typename parent_t::Flags;
  using typename parent_t::DomainHandlerType;
  using typename parent_t::topology_variant;
  using typename parent_t::eval_pts_variant;

  FormulaFunctionHandler(std::shared_ptr<FuncType> domain);

  virtual ~FormulaFunctionHandler() = default;

  void set_flags(const topology_variant &sdim,
                 const Flags &flag) override final;

  void fill_cache(const topology_variant &sdim,
                  ElementAccessor &elem,
                  const int s_id) const override final;

private:
  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const FuncType &func,
                        const self_t &func_handler,
                        ElementAccessor &elem,
                        const int s_id)
      :
      func_(func),
      func_handler_(func_handler),
      elem_(elem),
      s_id_(s_id)
    {}



    template<int sdim>
    void operator()(const Topology<sdim> &sub_elem)
    {
      using _D0 = function_element::template _D<0>;
      using _D1 = function_element::template _D<1>;
      using _D2 = function_element::template _D<2>;

      auto &local_cache = func_handler_.get_element_cache(elem_);
      auto &cache = local_cache.template get_sub_elem_cache<sdim>(s_id_);

      if (!cache.fill_none())
      {

        const auto &points = elem_.get_domain_element().template get_points<sdim>(s_id_);

        if (cache.template status_fill<_D0>())
        {
          auto &F = cache.template get_data<_D0>();
          func_.evaluate_0(points, F);
          F.set_status_filled(true);
        }
        if (cache.template status_fill<_D1>())
        {
          auto &DF = cache.template get_data<_D1>();
          func_.evaluate_1(points, DF);
          DF.set_status_filled(true);
        }
        if (cache.template status_fill<_D2>())
        {
          auto &D2F = cache.template get_data<_D2>();
          func_.evaluate_2(points, D2F);
          D2F.set_status_filled(true);
        }
//        if (cache.template status_fill<_Divergence>())
//          Assert(false,ExcNotImplemented());
      }

      cache.set_filled(true);
    }

    const FuncType &func_;
    const self_t     &func_handler_;
    ElementAccessor &elem_;
    const int s_id_;

  };

  friend struct FillCacheDispatcher;

};

IGA_NAMESPACE_CLOSE

#endif

