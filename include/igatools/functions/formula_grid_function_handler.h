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

#ifndef __FORMULA_GRID_FUNCTION_HANDLER_H_
#define __FORMULA_GRID_FUNCTION_HANDLER_H_

#include <igatools/functions/grid_function_handler.h>
#include <igatools/functions/formula_grid_function.h>

IGA_NAMESPACE_OPEN

/**
 *
 * @ingroup handlers
 */
template<int dim, int range>
class FormulaGridFunctionHandler :
  public GridFunctionHandler<dim, range>
{
private:
  using parent_t = GridFunctionHandler<dim, range>;
  using self_t = FormulaGridFunctionHandler<dim, range>;
protected:
  using typename parent_t::GridType;
public:
  using GridFunctionType =  const FormulaGridFunction<dim, range>;
  using typename parent_t::ElementAccessor;
  using typename parent_t::Flags;
  using typename parent_t::topology_variant;
  using typename parent_t::eval_pts_variant;

  FormulaGridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function);


  virtual ~FormulaGridFunctionHandler() = default;

  void set_flags(const topology_variant &sdim,
                 const Flags &flag) override final;

  void fill_cache(const topology_variant &sdim,
                  ElementAccessor &elem,
                  const int s_id) const override;

private:
  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    template <int order>
    using _D = typename ElementAccessor::template _D<order>;

    FillCacheDispatcher(const self_t &grid_function_handler,
                        ElementAccessor &elem,
                        const int s_id);


    template<int sdim>
    void operator()(const Topology<sdim> &sub_elem);

    const self_t     &grid_function_handler_;
    ElementAccessor &elem_;
    const int s_id_;
  };

//  friend struct FillCacheDispatcher;
};

IGA_NAMESPACE_CLOSE

#endif
