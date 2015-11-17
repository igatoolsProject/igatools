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

#ifndef __SUB_GRID_FUNCTION_HANDLER_H_
#define __SUB_GRID_FUNCTION_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/geometry/grid_function_handler.h>
#include <igatools/geometry/sub_grid_function.h>

IGA_NAMESPACE_OPEN




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
  SubGridFunctionHandler(const std::shared_ptr<GridFunctionType> &grid_function);


  virtual ~SubGridFunctionHandler() = default;




public:
  virtual void set_flags(const topology_variant &topology,
                         const Flags &flag) override;

  virtual void init_cache(GridFunctionElement<sdim,space_dim> &sub_grid_func_elem,
                          const eval_pts_variant &quad) const override;


  virtual void fill_cache(const topology_variant &topology,
                          GridFunctionElement<sdim,space_dim> &sub_grid_func_elem,
                          const int s_id) const override;



private:

  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const SubGridFunctionHandler<sdim,dim,space_dim> &sub_grid_func_handler,
                        SubGridFunctionElement<sdim,dim,space_dim> &sub_grid_func_elem,
                        const int s_id);

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

