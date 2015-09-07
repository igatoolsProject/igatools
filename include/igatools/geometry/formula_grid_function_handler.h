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

#ifndef __FORMULA_GRID_FUNCTION_HANDLER_H_
#define __FORMULA_GRID_FUNCTION_HANDLER_H_

#include <igatools/geometry/grid_function_handler.h>
#include <igatools/geometry/formula_grid_function.h>

IGA_NAMESPACE_OPEN

/**
 *
 */
template<int dim, int space_dim>
class FormulaGridFunctionHandler :
  public  GridFunctionHandler<dim, space_dim>
{
private:
  using parent_t = GridFunctionHandler<dim, space_dim>;
  using self_t = FormulaGridFunctionHandler<dim, space_dim>;
protected:
  using typename parent_t::GridType;
public:
  using GridFunctionType =  const FormulaGridFunction<dim,codim>;
  using typename parent_t::ConstElementAccessor;

  using typename parent_t::topology_variant;
  using typename parent_t::eval_pts_variant;

  FormulaGridFunctionHandler(std::shared_ptr<GridFunctionType> grid_function);


  virtual ~FormulaGridFunctionHandler() = default;


  void fill_cache(const topology_variant &sdim,
                  ConstElementAccessor &elem,
                  const int s_id) const override;

private:
  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const GridFunctionType &grid_function,
                        const self_t &grid_function_handler,
                        ConstElementAccessor &elem,
                        const int s_id)
      :
      grid_function_(grid_function),
      grid_function_handler_(grid_function_handler),
      elem_(elem),
      s_id_(s_id)
    {}


    template<int sdim>
    void operator()(const Topology<sdim> &sub_elem)
    {
      using _Point = typename ConstElementAccessor::_Point;
      using _Gradient = typename ConstElementAccessor::_Gradient;

      auto &local_cache = grid_function_handler_.get_element_cache(elem_);
      auto &cache = local_cache->template get_sub_elem_cache<sdim>(s_id_);

      if (!cache.fill_none())
      {
        const auto &grid_pts = elem_.get_grid_element().template get_points<sdim>(s_id_);
        if (cache.template status_fill<grid_function_element::_Point>())
        {
          grid_function_.evaluate_0(grid_pts, cache.template get_data<_Point>());
          cache.template set_status_filled<grid_function_element::_Point>(true);
        }

        if (cache.template status_fill<_Gradient>())
        {
          grid_function_.evaluate_1(grid_pts, cache.template get_data<_Gradient>());
          cache.template set_status_filled<_Gradient>(true);
        }
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

    const GridFunctionType &grid_function_;
    const self_t     &grid_function_handler_;
    ConstElementAccessor &elem_;
    const int s_id_;
  };

  friend struct FillCacheDispatcher;
};

IGA_NAMESPACE_CLOSE

#endif
