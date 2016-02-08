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

#ifndef __SUB_GRID_FUNCTION_HANDLER_H_
#define __SUB_GRID_FUNCTION_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/functions/grid_function_handler.h>
#include <igatools/functions/sub_grid_function.h>

IGA_NAMESPACE_OPEN




template<int sdim,int dim, int range>
class SubGridFunctionHandler
  : public GridFunctionHandler<sdim,range>
{
private:
  using parent_t = GridFunctionHandler<sdim,range>;
  using self_t = SubGridFunctionHandler<sdim,dim,range>;

public:
  using GridFunctionType = const SubGridFunction<sdim,dim,range>;

  using ElementAccessor = SubGridFunctionElement<sdim,dim,range>;
  using ElementIterator = GridIterator<ElementAccessor>;


  using Flags = grid_function_element::Flags;

protected:
  using FlagsArray = SafeSTLArray<Flags, sdim+1>;

  using topology_variant = TopologyVariants<sdim>;

  template<int k>
  using ConstQuad = const Quadrature<k>;
  using eval_pts_variant = SubElemPtrVariants<ConstQuad,sdim>;

private:
  SubGridFunctionHandler() = delete;

public:
  SubGridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function);


  virtual ~SubGridFunctionHandler() = default;




public:
  virtual void set_flags(const topology_variant &topology,
                         const Flags &flag) override;

  virtual void init_cache(GridFunctionElement<sdim,range> &sub_grid_func_elem,
                          const eval_pts_variant &quad) const override;


  virtual void fill_cache(const topology_variant &topology,
                          GridFunctionElement<sdim,range> &sub_grid_func_elem,
                          const int s_id) const override;



private:

  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const SubGridFunctionHandler<sdim,dim,range> &sub_grid_func_handler,
                        SubGridFunctionElement<sdim,dim,range> &sub_grid_func_elem,
                        const int s_id,
                        const int sup_grid_func_s_id);

    template<int k>
    void operator()(const Topology<k> &topology);


    const SubGridFunctionHandler<sdim,dim,range> &sub_grid_func_handler_;
    SubGridFunctionElement<sdim,dim,range> &sub_grid_func_elem_;
    const int s_id_;
    const int sup_grid_func_s_id_;
  };

private:

  std::unique_ptr<GridFunctionHandler<dim,range>> sup_grid_func_handler_;
  int sup_grid_func_s_id_;
};


IGA_NAMESPACE_CLOSE

#endif

