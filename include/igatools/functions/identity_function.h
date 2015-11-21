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

#ifndef __IDENTITY_FUNCTION_H_
#define __IDENTITY_FUNCTION_H_

#include <igatools/functions/function.h>
#include <igatools/geometry/grid.h>

IGA_NAMESPACE_OPEN
#if 0
template <int, int> class IdentityFunctionElementHandler;

/**
 * The identity function from R^dim to R^spacedim,
 * if dim < space_dim, the last space_dim-dim  coordinates are zero.
 *
 * @note this function is not inherited from formula function because
 * we want to optimize its computation
 *
 *
 * @ingroup serializable
 *
 * @author martinelli 2015
 * @author pauletti 2015
 */
template<int dim, int space_dim = dim>
class IdentityFunction :
  public Function<dim, 0, space_dim, 1>
{
private:
  using parent_t = Function<dim, 0, space_dim, 1>;
  using self_t = IdentityFunction<dim,space_dim>;

protected:
  using GridType = Grid<dim>;

public:
  using typename parent_t::ElementAccessor;
  using ElementHandler = IdentityFunctionElementHandler<dim, space_dim>;


private:
  IdentityFunction();


public:
  IdentityFunction(std::shared_ptr<const GridType> grid);
  virtual ~IdentityFunction() = default;

  static std::shared_ptr<parent_t>
  create(std::shared_ptr<const GridType> grid)
  {
    return std::make_shared<self_t>(grid);
  }

  static std::shared_ptr<const parent_t>
  const_create(std::shared_ptr<const GridType> grid)
  {
    return std::make_shared<const self_t>(grid);
  }

  std::shared_ptr<typename parent_t::ElementHandler>
  create_cache_handler() const override final;

private:
  std::shared_ptr<const GridType> grid_;

#ifdef MESH_REFINEMENT

  void create_connection_for_insert_knots(std::shared_ptr<self_t> &ig_function);

  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid);

#endif // MESH_REFINEMENT

#ifdef SERIALIZATION
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class boost::serialization::access;

  template<class Archive>
  void
  serialize(Archive &ar, const unsigned int version)
  {
    ar &boost::serialization::make_nvp("IdentityFunction_base_t",
                                       boost::serialization::base_object<parent_t>(*this));
  }
  ///@}
#endif // SERIALIZATION

};
#endif
IGA_NAMESPACE_CLOSE

#endif
