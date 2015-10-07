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

#ifndef __IG_GRID_FUNCTION_HANDLER_H_
#define __IG_GRID_FUNCTION_HANDLER_H_

#include <igatools/geometry/grid_function_handler.h>
#include <igatools/functions/ig_grid_function.h>
#include <igatools/basis_functions/space_element_handler.h>


IGA_NAMESPACE_OPEN



/**
 * @ingroup handlers
 */
template<int dim, int space_dim>
class IgGridFunctionHandler :
  public GridFunctionHandler<dim, space_dim>
{
private:
  using parent_t = GridFunctionHandler<dim, space_dim>;
  using self_t = IgGridFunctionHandler<dim, space_dim>;
protected:
  using typename parent_t::GridType;
public:
  using GridFunctionType =  const IgGridFunction<dim, space_dim>;
  using typename parent_t::ConstElementAccessor;
  using typename parent_t::Flags;
  using typename parent_t::topology_variant;
  using typename parent_t::eval_pts_variant;

  IgGridFunctionHandler(const std::shared_ptr<GridFunctionType> &ig_grid_function);


  virtual ~IgGridFunctionHandler() = default;

  void set_flags(const topology_variant &sdim,
                 const Flags &flag) override final;

  void fill_cache(const topology_variant &sdim,
                  ConstElementAccessor &elem,
                  const int s_id) const override final;

private:
  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    template <int order>
    using _D = typename ConstElementAccessor::template _D<order>;

    FillCacheDispatcher(const GridFunctionType &ig_grid_function,
                        const self_t &ig_grid_function_handler,
                        const GridHandler<dim> &grid_handler,
                        ConstElementAccessor &ig_grid_function_elem,
                        const int s_id)
      :
      ig_grid_function_(ig_grid_function),
      ig_grid_function_handler_(ig_grid_function_handler),
      grid_handler_(grid_handler),
      ig_grid_function_elem_(ig_grid_function_elem),
      s_id_(s_id)
    {}


    template<int sdim>
    void operator()(const Topology<sdim> &sub_elem);


    const GridFunctionType &ig_grid_function_;
    const self_t &ig_grid_function_handler_;
    const GridHandler<dim> &grid_handler_;
    ConstElementAccessor &ig_grid_function_elem_;
    const int s_id_;
  };

//  friend struct FillCacheDispatcher;

  using IgSpaceHandler = SpaceElementHandler<dim,0,space_dim,1,Transformation::h_grad>;
  std::unique_ptr<IgSpaceHandler> ig_space_handler_;

};

IGA_NAMESPACE_CLOSE

#endif
