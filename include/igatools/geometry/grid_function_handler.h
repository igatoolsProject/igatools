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

protected:
  std::unique_ptr<GridHandler> grid_handler_;

private:
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

  using ElementAccessor = SubGridFunctionElement<sdim,dim,space_dim>;
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
  SubGridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function)
    :
    parent_t(grid_function),
    sup_grid_func_handler_(grid_function->get_sup_grid_function()->create_cache_handler())
  {}


  virtual ~SubGridFunctionHandler() = default;




public:
  virtual void set_flags(const topology_variant &topology,
                         const Flags &flag) override
  {
    parent_t::set_flags(topology,flag);

    sup_grid_func_handler_->set_flags(topology,flag);
  }

  virtual void init_cache(GridFunctionElement<sdim,space_dim> &sub_grid_func_elem,
                          const eval_pts_variant &quad) const override
  {
    parent_t::init_cache(sub_grid_func_elem,quad);

    auto &as_sub_grid_func_elem =
      dynamic_cast<SubGridFunctionElement<sdim,dim,space_dim> &>(sub_grid_func_elem);

    this->sup_grid_func_handler_->init_cache(as_sub_grid_func_elem.get_sup_grid_function_element(),quad);
  }


  virtual void fill_cache(const topology_variant &topology,
                          GridFunctionElement<sdim,space_dim> &sub_grid_func_elem,
                          const int s_id) const override
  {
    this->grid_handler_->fill_cache(topology, sub_grid_func_elem.get_grid_element(), s_id);

    auto &as_sub_grid_func_elem =
      dynamic_cast<SubGridFunctionElement<sdim,dim,space_dim> &>(sub_grid_func_elem);

    this->sup_grid_func_handler_->fill_cache(topology,as_sub_grid_func_elem.get_sup_grid_function_element(),s_id);

    auto fill_dispatcher = FillCacheDispatcher(
                             *this,
                             as_sub_grid_func_elem,
                             s_id);
    boost::apply_visitor(fill_dispatcher, topology);
  }



private:

  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const SubGridFunctionHandler<sdim,dim,space_dim> &sub_grid_func_handler,
                        SubGridFunctionElement<sdim,dim,space_dim> &sub_grid_func_elem,
                        const int s_id)
      :
      sub_grid_func_handler_(sub_grid_func_handler),
      sub_grid_func_elem_(sub_grid_func_elem),
      s_id_(s_id)
    {}

    template<int k>
    void operator()(const Topology<k> &topology)
    {
      static_assert(k >= 0 && k <= sdim,"Invalid topological dimension.");


      const auto &sub_unit_elem = UnitElement<dim>::template get_elem<k>(s_id_);
      const auto &sub_elem_active_dirs = sub_unit_elem.active_directions;


      auto &sub_grid_func_local_cache = sub_grid_func_handler_.get_element_cache(sub_grid_func_elem_);
      auto &sub_grid_func_cache = sub_grid_func_local_cache.template get_sub_elem_cache<k>(s_id_);

      if (!sub_grid_func_cache.fill_none())
      {
        const auto &sup_grid_func_elem = sub_grid_func_elem_.get_sup_grid_function_element();

        using _D0 = typename grid_function_element::template _D<0>;
        if (sub_grid_func_cache.template status_fill<_D0>())
        {
          const auto &sup_grid_func_D0 = sup_grid_func_elem.template get_values_from_cache<_D0,k>(s_id_);
          auto &D0 = sub_grid_func_cache.template get_data<_D0>();

          D0 = sup_grid_func_D0;
          D0.set_status_filled(true);
        }

        using _D1 = typename grid_function_element::template _D<1>;
        if (sub_grid_func_cache.template status_fill<_D1>())
        {
          const auto &sup_grid_func_D1 = sup_grid_func_elem.template get_values_from_cache<_D1,k>(s_id_);
          auto &D1 = sub_grid_func_cache.template get_data<_D1>();

          const int n_pts = sup_grid_func_D1.get_num_points();
          for (int pt = 0 ; pt < n_pts ; ++pt)
          {
            auto &D1_pt = D1[pt];
            const auto sup_grid_func_D1_pt = sup_grid_func_D1[pt];

            int dir = 0;
            for (const int active_dir : sub_elem_active_dirs)
            {
              D1_pt[dir] = sup_grid_func_D1_pt[active_dir];
              ++dir;
            }
          } // end loop pt

          D1.set_status_filled(true);
        }

        using _D2 = typename grid_function_element::template _D<2>;
        if (sub_grid_func_cache.template status_fill<_D2>())
        {
          const auto &sup_grid_func_D2 = sup_grid_func_elem.template get_values_from_cache<_D2,k>(s_id_);
          auto &D2 = sub_grid_func_cache.template get_data<_D2>();

          const int n_pts = sup_grid_func_D2.get_num_points();
          for (int pt = 0 ; pt < n_pts ; ++pt)
          {
            auto &D2_pt = D2[pt];
            const auto sup_grid_func_D2_pt = sup_grid_func_D2[pt];

            int dir_i = 0;
            for (const int active_dir_i : sub_elem_active_dirs)
            {
              auto &D2_pt_i = D2_pt[dir_i];
              const auto &sup_grid_func_D2_pt_i = sup_grid_func_D2_pt[active_dir_i];

              int dir_j = 0;
              for (const int active_dir_j : sub_elem_active_dirs)
              {
                D2_pt_i[dir_j] = sup_grid_func_D2_pt_i[active_dir_j];
                ++dir_j;
              }
              ++dir_i;
            }
          } // end loop pt

          D2.set_status_filled(true);
        }
//        if (cache.template status_fill<_Divergence>())
//          Assert(false,ExcNotImplemented());
      }

      sub_grid_func_cache.set_filled(true);
    }




    const SubGridFunctionHandler<sdim,dim,space_dim> &sub_grid_func_handler_;
    SubGridFunctionElement<sdim,dim,space_dim> &sub_grid_func_elem_;
    const int s_id_;
  };

private:

  std::unique_ptr<GridFunctionHandler<dim,space_dim>> sup_grid_func_handler_;
};


IGA_NAMESPACE_CLOSE

#endif

