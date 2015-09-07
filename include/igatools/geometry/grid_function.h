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

#ifndef __GRID_FUNCTION_H_
#define __GRID_FUNCTION_H_

#include <igatools/base/config.h>
#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_handler.h>

IGA_NAMESPACE_OPEN

template <int, int, class> class GridFunctionElementBase;
template <int, int> class GridFunctionElement;
template <int, int> class ConstGridFunctionElement;
template <int, int> class GridFunctionHandler;

/**
 *
 */
template<int dim_, int space_dim_>
class GridFunction :
  public std::enable_shared_from_this<GridFunction<dim_,space_dim_> >
{
private:
  using self_t = GridFunction<dim_, space_dim_>;

public:
  static const int space_dim = space_dim_;
  static const int dim = dim_;

  using GridType = const Grid<dim_>;

  using ElementAccessor = GridFunctionElement<dim_, codim_>;
  using ElementIterator = GridIterator<ElementAccessor>;
  using ConstElementAccessor = ConstGridFunctionElement<dim_, codim_>;
  using ElementConstIterator = GridIterator<ConstElementAccessor>;

  using ElementHandler = GridFunctionHandler<dim_, codim_>;

  using List = typename GridType::List;
  using ListIt = typename GridType::ListIt;

public:
  using GridPoint = typename GridType::Point;
  using Value = Values<dim, space_dim, 1>;
  template <int order>
  using Derivative = Derivatives<dim, space_dim, 1, order>;

  using Gradient = Derivative<1>;
  ///@}

private:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  GridFunction() = default;

protected:
  GridFunction(std::shared_ptr<GridType> grid);

public:
  ~GridFunction();

  static std::shared_ptr<self_t>
  create(std::shared_ptr<GridType> grid)
  {
    return std::shared_ptr<self_t>(new self_t(grid));
  }


  static std::shared_ptr<const self_t>
  const_create(std::shared_ptr<GridType> grid)
  {
    return create(grid);
  }

  std::shared_ptr<GridType> get_grid() const;

public:
  virtual std::shared_ptr<ElementHandler>
  create_cache_handler() const;

  std::unique_ptr<ConstElementAccessor>
  create_element(const ListIt &index, const PropId &prop) const;

  std::unique_ptr<ElementAccessor>
  create_element(const ListIt &index, const PropId &prop);

  ///@name Iterating of grid elements
  ///@{
  /**
   * This function returns a element iterator to the first element of the patch.
   */
  ElementIterator begin(const PropId &prop = ElementProperties::active);

  /**
   * This function returns a element iterator to one-pass the end of patch.
   */
  ElementIterator end(const PropId &prop = ElementProperties::active);

  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  ElementConstIterator begin(const PropId &prop = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  ElementConstIterator end(const PropId &prop = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  ElementConstIterator cbegin(const PropId &prop = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  ElementConstIterator cend(const PropId &prop = ElementProperties::active) const;
  ///@}


  void print_info(LogStream &out) const
  {
    AssertThrow(false,ExcNotImplemented());
  }

private:
  std::shared_ptr<GridType> grid_;

  friend class GridFunctionElementBase<dim_, codim_, GridFunction<dim_, codim_>>;
  friend class GridFunctionElementBase<dim_, codim_, const GridFunction<dim_, codim_>>;
  friend class GridFunctionElement<dim_, codim_>;
  friend class ConstGridFunctionElement<dim_, codim_>;

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
    ar.template register_type<IgFunction<dim_,0,dim_+codim_,1> >();
    ar &boost::serialization::make_nvp("F_",F_);
    ar &boost::serialization::make_nvp("flags_",flags_);
  }
  ///@}
#endif
};

IGA_NAMESPACE_CLOSE

#endif

