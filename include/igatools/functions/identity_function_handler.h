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

#ifndef __IDENTITY_FUNCTION_HANDLER_H_
#define __IDENTITY_FUNCTION_HANDLER_H_

#include <igatools/functions/function.h>
#include <igatools/functions/function_handler.h>
#include <igatools/geometry/grid_handler.h>

IGA_NAMESPACE_OPEN
#if 0
/**
 * create a rectangular "identity" tensor
 */
template<int dim,int space_dim>
auto
create_id_tensor()
{
  typename Function<dim, 0, space_dim, 1>::Gradient res;
  for (int i=0; i<dim; ++i)
    res[i][i] = 1.;
  return res;
}


/**
 * The identity function from R^dim to R^spacedim,
 * if dim < space_dim, the last space_dim-dim  coordinates are zero.
 *
 * @note this function is not inherited from formula function because
 * we want to optimize its computation
 *
 *
 *
 * @author martinelli 2015
 * @author pauletti 2015
 */
template<int dim, int space_dim = dim>
class IdentityFunctionHandler :
  public FunctionHandler<dim, 0, space_dim, 1>
{
private:
  using parent_t = FunctionHandler<dim, 0, space_dim, 1>;
  using self_t = IdentityFunctionHandler<dim,space_dim>;

protected:
  using GridType = Grid<dim>;
  using GridHandlerType = typename GridType::Handler;

  using FuncType = Function<dim, 0, space_dim, 1>;

public:
  using typename parent_t::Flags;
  using typename parent_t::ConstElementAccessor;
  using typename parent_t::eval_pts_variant;
  using typename parent_t::topology_variant;
#if 0
  using typename parent_t::ElementAccessor;
  using typename parent_t::_Value;
  using typename parent_t::_Gradient;
// using typename parent_t::_Hessian;

  using typename parent_t::topology_variant;
  //using typename parent_t::Point;
  using typename parent_t::Value;
  using typename parent_t::Gradient;
  using typename parent_t::Hessian;
#endif
  //  using typename parent_t::ElementIterator;
  // using typename parent_t::ElementAccessor;

//    template <int order>
//    using Derivative = typename parent_t::template Derivative<order>;
  // IdentityFunctionHandler(std::shared_ptr<GridType> grid);

private:
  /**
   * Default constructor. Not allowed to be used.
   */
  IdentityFunctionHandler() = delete;

public:
  IdentityFunctionHandler(std::shared_ptr<const GridType> grid);

  virtual ~IdentityFunctionHandler() = default;

  static std::shared_ptr<parent_t>
  create(std::shared_ptr<const GridType> grid)
  {
    return std::make_shared<self_t>(grid);
  }



  void set_flags(const topology_variant &sdim,
                 const Flags &flag) override final;

  void init_cache(ConstElementAccessor &elem,
                  const eval_pts_variant &quad) const override final;

  void fill_cache(const topology_variant &sdim,
                  ConstElementAccessor &elem,
                  const int s_id) const override final;




private:
  std::shared_ptr<GridHandlerType> grid_handler_;
#if 0
  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const int s_id, const self_t &func, ElementAccessor &elem)
      :
      s_id_(s_id),
      funct_(func),
      elem_(elem)
    {}

    template<int sdim>
    void operator()(const Topology<sdim> &sub_elem)
    {
      auto &local_cache = funct_.get_cache(elem_);
      auto &cache = local_cache->template get_sub_elem_cache<sdim>(s_id_);

      if (!cache.fill_none())
      {
        if (cache.template status_fill<_Value>())
        {
          const auto points =
            elem_.get_domain_element()->
            get_grid_element()->template get_points<sdim>(s_id_);
          if (cache.template status_fill<_Value>())
          {
            auto &values = cache.template get_data<_Value>();
            values = points;
            cache.template set_status_filled<_Value>(true);
          }
        }
        if (cache.template status_fill<_Gradient>())
        {
          // TODO (pauletti, Apr 17, 2015): this can be static const
          const auto identity = create_id_tensor<dim,space_dim>();
          cache.template get_data<_Gradient>().fill(identity);

          cache.template set_status_filled<_Gradient>(true);
        }
//        if (cache.template status_fill<_Hessian>())
//        {
//          // TODO (pauletti, Apr 17, 2015): this can be static const
//          Hessian zero;
//          cache.template get_data<_Hessian>().fill(zero);
//
//          cache.template set_status_filled<_Hessian>(true);
//        }
      }
      cache.set_filled(true);
    }

    const int s_id_;
    const self_t &funct_;
    ElementAccessor &elem_;
  };

  friend struct FillCacheDispatcher;
#endif

#ifdef MESH_REFINEMENT

  void create_connection_for_insert_knots(std::shared_ptr<self_t> &ig_function);

  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid);

#endif // MESH_REFINEMENT

#if 0
#ifdef SERIALIZATION
  /**
   * @name FunctionHandlers needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class boost::serialization::access;

  template<class Archive>
  void
  serialize(Archive &ar, const unsigned int version)
  {
    ar &boost::serialization::make_nvp("IdentityFunctionHandler_base_t",
                                       boost::serialization::base_object<parent_t>(*this));
  }
  ///@}
#endif // SERIALIZATION
#endif
};

#endif
IGA_NAMESPACE_CLOSE

#endif
