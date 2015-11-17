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


  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;
  using Flags = grid_function_element::Flags;

protected:
  using FlagsArray = SafeSTLArray<Flags, dim+1>;

  using topology_variant = TopologyVariants<dim_>;

  template<int k>
  using ConstQuad = const Quadrature<k>;
  using eval_pts_variant = SubElemPtrVariants<ConstQuad,dim_>;

private:
  GridFunctionHandler() = delete;

public:
  GridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function);


  virtual ~GridFunctionHandler() = default;


  std::shared_ptr<GridFunctionType> get_grid_function() const;


public:
  virtual void set_flags(const topology_variant &sdim,
                         const Flags &flag);

  template <int sdim>
  void set_flags(const Flags &flag)
  {
    this->set_flags(Topology<sdim>(), flag);
  }

  void set_element_flags(const Flags &flag);

  virtual void init_cache(ElementAccessor &elem,
                          const eval_pts_variant &quad) const;

  void init_cache(ElementIterator &elem,
                  const eval_pts_variant &quad) const;

  void init_element_cache(ElementAccessor &elem,
                          const std::shared_ptr<const Quadrature<dim_>> &quad) const;

  void init_element_cache(ElementIterator &elem,
                          const std::shared_ptr<const Quadrature<dim_>> &quad) const;

  virtual void fill_cache(const topology_variant &sdim,
                          ElementAccessor &elem,
                          const int s_id) const = 0;

  void fill_cache(const topology_variant &sdim,
                  ElementIterator &elem,
                  const int s_id) const;

  template <int sdim>
  void fill_cache(ElementIterator &elem,
                  const int s_id) const
  {
    this->fill_cache(Topology<sdim>(), elem, s_id);
  }

  template <int sdim>
  void fill_cache(ElementAccessor &elem,
                  const int s_id) const
  {
    this->fill_cache(Topology<sdim>(), elem, s_id);
  }

  void fill_element_cache(ElementAccessor &elem) const;

  void fill_element_cache(ElementIterator &elem) const;

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
    SetFlagsDispatcher(const Flags flag, FlagsArray &flags)
      :
      flag_(flag),
      flags_(flags)
    {}

    template<int sdim>
    void operator()(const Topology<sdim> &)
    {
      flags_[sdim] |= flag_;
    }

    const Flags flag_;
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



template<int sdim,int dim, int space_dim>
class SubGridFunctionHandler
  : public GridFunctionHandler<sdim,space_dim>
{
private:
  using parent_t = GridFunctionHandler<sdim,space_dim>;
  using self_t = SubGridFunctionHandler<sdim,dim,space_dim>;

public:
  using GridFunctionType = const SubGridFunction<sdim,dim,space_dim>;
//  using GridType = const Grid<dim_>;
//  using GridHandler = typename GridType::ElementHandler;

  using ElementAccessor = SubGridFunctionElement<sdim,dim,space_dim>;
  using ElementIterator = GridIterator<ElementAccessor>;


//  using List = typename GridType::List;
//  using ListIt = typename GridType::ListIt;
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
  SubGridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function)
    :
    parent_t(grid_function),
    sup_grid_func_handler_(grid_function->get_sup_grid_function()->create_cache_handler())
  {}


  virtual ~SubGridFunctionHandler() = default;




public:
  virtual void set_flags(const topology_variant &topology,
                         const Flags &flag)
  {
    parent_t::set_flags(topology,flag);

    auto set_flags_dispatcher = SetFlagsDispatcher(flag,*sup_grid_func_handler_);
    boost::apply_visitor(set_flags_dispatcher,topology);
  }

  virtual void init_cache(ElementAccessor &sub_grid_func_elem,
                          const eval_pts_variant &quad) const
  {
    parent_t::init_cache(sub_grid_func_elem,quad);

    auto init_dispatcher = InitCacheDispatcher(
                             *(this->sup_grid_func_handler_),
                             sub_grid_func_elem.get_sup_grid_function_element());
    boost::apply_visitor(init_dispatcher, quad);
  }


  virtual void fill_cache(const topology_variant &topology,
                          GridFunctionElement<sdim,space_dim> &sub_grid_func_elem,
                          const int s_id) const override
  {
    auto &as_sub_grid_func_elem =
      dynamic_cast<SubGridFunctionElement<sdim,dim,space_dim> &>(sub_grid_func_elem);

    auto fill_dispatcher = FillCacheDispatcher(
                             *this,
                             *(this->sup_grid_func_handler_),
                             as_sub_grid_func_elem,
                             as_sub_grid_func_elem.get_sup_grid_function_element(),
                             s_id);
    boost::apply_visitor(fill_dispatcher, topology);
  }



private:
  /**
   * Alternative to
   * template <int k> set_flags()
   */
  struct SetFlagsDispatcher : boost::static_visitor<void>
  {
    SetFlagsDispatcher(
      const Flags flag,
      GridFunctionHandler<dim,space_dim> &sup_handler)
      :
      flag_(flag),
      sup_handler_(sup_handler)
    {}

    template<int k>
    void operator()(const Topology<k> &topology)
    {
      Assert(k >= 0 && k <= sdim,ExcMessage("Invalid topology dimension."));
      sup_handler_.template set_flags<k>(flag_);
    }

    const Flags flag_;
    GridFunctionHandler<dim,space_dim> &sup_handler_;
  };


  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const GridFunctionHandler<dim,space_dim> &sup_grid_func_handler,
                        GridFunctionElement<dim,space_dim> &sup_grid_func_elem)
      :
      sup_grid_func_handler_(sup_grid_func_handler),
      sup_grid_func_elem_(sup_grid_func_elem)
    {}


    template<int k>
    void operator()(const std::shared_ptr<const Quadrature<k>> &quad)
    {
      Assert(k >= 0 && k <= sdim,ExcMessage("Invalid topology dimension for the Quadrature scheme."));
      sup_grid_func_handler_.init_cache(sup_grid_func_elem_,quad);
    }

    const GridFunctionHandler<dim,space_dim> &sup_grid_func_handler_;
    GridFunctionElement<dim,space_dim> &sup_grid_func_elem_;
//    const FlagsArray &flags_;
  };



  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const SubGridFunctionHandler<sdim,dim,space_dim> &sub_grid_func_handler,
                        const GridFunctionHandler<dim,space_dim> &sup_grid_func_handler,
                        SubGridFunctionElement<sdim,dim,space_dim> &sub_grid_func_elem,
                        GridFunctionElement<dim,space_dim> &sup_grid_func_elem,
                        const int s_id)
      :
      sub_grid_func_handler_(sub_grid_func_handler),
      sup_grid_func_handler_(sup_grid_func_handler),
      sub_grid_func_elem_(sub_grid_func_elem),
      sup_grid_func_elem_(sup_grid_func_elem),
      s_id_(s_id)
    {}

    template<int k>
    void operator()(const Topology<k> &topology)
    {
      Assert(k >= 0 && k <= sdim,ExcMessage("Invalid topology dimension."));
      sup_grid_func_handler_.template fill_cache<k>(sup_grid_func_elem_,s_id_);



      auto &sub_grid_func_local_cache = sub_grid_func_handler_.get_element_cache(sub_grid_func_elem_);
      auto &sub_grid_func_cache = sub_grid_func_local_cache.template get_sub_elem_cache<k>(s_id_);

      if (!sub_grid_func_cache.fill_none())
      {
//        const auto &grid_pts = elem_.get_grid_element().template get_points<sdim>(s_id_);
        using _D0 = typename grid_function_element::template _D<0>;
        if (sub_grid_func_cache.template status_fill<_D0>())
        {
          Assert(false,ExcNotImplemented());

//          auto &F = cache.template get_data<_D<0>>();
//          formula_grid_function.evaluate_0(grid_pts, F);
//          F.set_status_filled(true);
        }

        using _D1 = typename grid_function_element::template _D<1>;
        if (sub_grid_func_cache.template status_fill<_D1>())
        {
          Assert(false,ExcNotImplemented());
//          auto &DF = cache.template get_data<_D<1>>();
//          formula_grid_function.evaluate_1(grid_pts, DF);
//          DF.set_status_filled(true);
        }

        using _D2 = typename grid_function_element::template _D<2>;
        if (sub_grid_func_cache.template status_fill<_D2>())
        {
          Assert(false,ExcNotImplemented());
//          auto &D2F = cache.template get_data<_D<2>>();
//          formula_grid_function.evaluate_2(grid_pts, D2F);
//          D2F.set_status_filled(true);
        }
//        if (cache.template status_fill<_Divergence>())
//          Assert(false,ExcNotImplemented());
      }

      sub_grid_func_cache.set_filled(true);
    }




    const SubGridFunctionHandler<sdim,dim,space_dim> &sub_grid_func_handler_;
    const GridFunctionHandler<dim,space_dim> &sup_grid_func_handler_;
    SubGridFunctionElement<sdim,dim,space_dim> &sub_grid_func_elem_;
    GridFunctionElement<dim,space_dim> &sup_grid_func_elem_;
    const int s_id_;
//    const FlagsArray &flags_;
  };

private:

  std::unique_ptr<GridFunctionHandler<dim,space_dim>> sup_grid_func_handler_;
//  std::unique_ptr<GridHandler> grid_handler_;

//  FlagsArray flags_;

//  friend ElementAccessor;
};


IGA_NAMESPACE_CLOSE

#endif

