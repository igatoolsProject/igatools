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

#ifndef __GRID_FUNCTION_HANDLER_H_
#define __GRID_FUNCTION_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_function.h>
#include <igatools/geometry/grid_handler.h>

IGA_NAMESPACE_OPEN

template <int, int> class GridFunctionElement;

/**
 *
 * @ingroup handlers
 */
template<int dim_, int space_dim_>
class GridFunctionHandler
//    :
//  public std::enable_shared_from_this<GridFunctionHandler<dim_,space_dim_> >
{
private:
  using self_t = GridFunctionHandler<dim_, space_dim_>;

public:
  static const int space_dim = space_dim_;
  static const int dim = dim_;

  using GridFunctionType = const GridFunction<dim_, space_dim_>;
  using GridType = const Grid<dim_>;
  using GridHandler = typename GridType::ElementHandler;

  using ElementAccessor = GridFunctionElement<dim_, space_dim_>;
  using ElementIterator = GridIterator<ElementAccessor>;
//  using ConstElementAccessor = GridFunctionElement<dim_, space_dim_>;
  using ElementConstIterator = GridIterator<ElementAccessor>;

  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;
  using Flags = grid_function_element::Flags;
  using CacheFlags = grid_function_element::CacheFlags;

protected:
  using FlagsArray = SafeSTLArray<CacheFlags, dim+1>;

  using topology_variant = TopologyVariants<dim_>;

  template<int k>
  using ConstQuad = const Quadrature<k>;
  using eval_pts_variant = SubElemPtrVariants<ConstQuad,dim_>;

private:
  GridFunctionHandler() = delete;

public:
  GridFunctionHandler(std::shared_ptr<GridFunctionType> grid_function);


  virtual ~GridFunctionHandler();

#if 0
  static std::shared_ptr<self_t>
  create(std::shared_ptr<GridFunctionType> grid_function)
  {
    return std::shared_ptr<self_t>(new self_t(grid_function));
  }


  static std::shared_ptr<const self_t>
  const_create(std::shared_ptr<GridFunctionType> grid_function)
  {
    return create(grid_function);
  }
#endif

  std::shared_ptr<GridFunctionType> get_grid_function() const;


public:
  virtual void set_flags(const topology_variant &sdim,
                         const Flags &flag);

  template <int sdim>
  void set_flags(const Flags &flag)
  {
    this->set_flags(Topology<sdim>(), flag);
  }

  virtual void init_cache(ElementAccessor &elem,
                          const eval_pts_variant &quad) const;

  void init_cache(ElementConstIterator &elem,
                  const eval_pts_variant &quad) const
  {
    this->init_cache(*elem, quad);
  }

  virtual void fill_cache(const topology_variant &sdim,
                          ElementAccessor &elem,
                          const int s_id) const;

  void fill_cache(const topology_variant &sdim,
                  ElementConstIterator &elem,
                  const int s_id) const
  {
    this->fill_cache(sdim, *elem, s_id);
  }

  template <int sdim>
  void fill_cache(ElementConstIterator &elem,
                  const int s_id)
  {
    this->fill_cache(Topology<sdim>(), elem, s_id);
  }

  template <int sdim>
  void fill_cache(ElementAccessor &elem,
                  const int s_id)
  {
    this->fill_cache(Topology<sdim>(), elem, s_id);
  }

  //protected:
public:
  const GridHandler &
  get_grid_handler() const;

  GridHandler &
  get_grid_handler();

protected:
  typename ElementAccessor::CacheType &
  get_element_cache(ElementAccessor &elem) const
  {
    return  elem.local_cache_;
  }
//*/

private:
  /**
   * Alternative to
   * template <int sdim> set_flags()
   */
  struct SetFlagsDispatcher : boost::static_visitor<void>
  {
    SetFlagsDispatcher(const CacheFlags flag, FlagsArray &flags)
      :
      flag_(flag),
      flags_(flags)
    {}

    template<int sdim>
    void operator()(const Topology<sdim> &)
    {
      flags_[sdim] |= flag_;
    }

    const CacheFlags flag_;
    FlagsArray &flags_;
  };



  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const self_t &grid_function_handler,
                        ElementAccessor &elem,
                        const FlagsArray &flags)
      :
      grid_function_handler_(grid_function_handler),
      elem_(elem),
      flags_(flags)
    {}


    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
    {
      auto &cache = grid_function_handler_.get_element_cache(elem_);

      const auto n_points = elem_.get_grid_element().template get_quad<sdim>()
                            ->get_num_points();
      for (auto &s_id: UnitElement<dim_>::template elems_ids<sdim>())
      {
        auto &s_cache = cache.template get_sub_elem_cache<sdim>(s_id);
        s_cache.resize(flags_[sdim], n_points);
      }
    }

    const self_t &grid_function_handler_;
    ElementAccessor &elem_;
    const FlagsArray &flags_;
  };




private:
  std::shared_ptr<GridFunctionType> grid_function_;

  std::unique_ptr<GridHandler> grid_handler_;

  FlagsArray flags_;

//  friend ElementAccessor;
};

IGA_NAMESPACE_CLOSE

#endif

