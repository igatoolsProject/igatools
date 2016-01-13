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

#include <igatools/functions/grid_function_handler.h>
#include <igatools/functions/ig_grid_function.h>
#include <igatools/basis_functions/space_element_handler.h>


IGA_NAMESPACE_OPEN



/**
 * @ingroup handlers
 */
template<int dim, int range>
class IgGridFunctionHandler :
  public GridFunctionHandler<dim, range>
{
private:
  using parent_t = GridFunctionHandler<dim, range>;
  using self_t = IgGridFunctionHandler<dim, range>;
protected:
  using typename parent_t::GridType;
public:
  using GridFunctionType = const IgGridFunction<dim, range>;
  using typename parent_t::ElementAccessor;
  using typename parent_t::Flags;
  using typename parent_t::topology_variant;
  using typename parent_t::eval_pts_variant;

  IgGridFunctionHandler(const std::shared_ptr<GridFunctionType> &ig_grid_function);


  virtual ~IgGridFunctionHandler() = default;

  void set_flags(const topology_variant &sdim,
                 const Flags &flag) override final;

  void fill_cache(const topology_variant &sdim,
                  ElementAccessor &elem,
                  const int s_id) const override final;


  std::shared_ptr<GridFunctionType>
  get_ig_grid_function() const;

private:
  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    template <int order>
    using _D = typename ElementAccessor::template _D<order>;

    FillCacheDispatcher(const self_t &ig_grid_function_handler,
                        ElementAccessor &ig_grid_function_elem,
                        const int s_id)
      :
      ig_grid_function_handler_(ig_grid_function_handler),
      ig_grid_function_elem_(ig_grid_function_elem),
      s_id_(s_id)
    {}


    template<int sdim>
    void operator()(const Topology<sdim> &sub_elem);


    const self_t &ig_grid_function_handler_;
    ElementAccessor &ig_grid_function_elem_;
    const int s_id_;
  };

//  friend struct FillCacheDispatcher;

  using IgBasisHandler = SpaceElementHandler<dim,0,range,1>;
  std::unique_ptr<IgBasisHandler> ig_basis_handler_;

  std::shared_ptr<GridFunctionType> ig_grid_function_;

};

IGA_NAMESPACE_CLOSE

#endif
