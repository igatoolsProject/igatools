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

#ifndef __DOMAIN_H_
#define __DOMAIN_H_

#include <igatools/base/config.h>
#include <igatools/geometry/grid_function.h>
#include <igatools/geometry/grid_function_handler.h>

IGA_NAMESPACE_OPEN

template <int, int, class> class DomainElementBase;
template <int, int> class DomainElement;
template <int, int> class ConstDomainElement;
template <int, int> class DomainHandler;

/**
 * @brief The mapping is a deformation \f$ F : \hat\Omega \to \Omega\f$
 * which maps the reference domain \f$\hat\Omega \in \mathbb{R}^{dim}\f$ to the
 * physical domain \f$\Omega \in \mathbb{R}^{dim+codim}\f$.
 *
 * Domain is the physical domain.
 *
 *
 *
 * @ingroup containers
 * @ingroup serializable
 *
 * @author pauletti 2014, 2015
 * @author M. Martinelli, 2015
 */
template<int dim_, int codim_ = 0>
class Domain :
  public std::enable_shared_from_this<Domain<dim_,codim_> >
{
private:
  using self_t = Domain<dim_, codim_>;

public:
  static const int space_dim = dim_ + codim_;
  static const int dim = dim_;

  using GridFuncType = const GridFunction<dim, space_dim>;

  using ElementAccessor = DomainElement<dim_, codim_>;
  using ElementIterator = GridIterator<ElementAccessor>;
  using ConstElementAccessor = ConstDomainElement<dim_, codim_>;
  using ElementConstIterator = GridIterator<ConstElementAccessor>;

  using ElementHandler = DomainHandler<dim_, codim_>;

  using List = typename GridFuncType::List;
  using ListIt = typename GridFuncType::ListIt;

public:
  //using GridPoint = typename GridFuncType::Point;
  using Point =  Points<space_dim>;
  template <int order>
  using Derivative = Derivatives<dim, space_dim, 1, order>;


//  /** Type for the return of the function.*/
//  using Value = Values<dim, space_dim, 1>;
//
//  /**
//   * Type for the derivative of the function.
//   */

//
//  /**
//   * Type for the gradient of the function.
//   */
  using Gradient = Derivative<1>;
//
//  /**
//   * Type for the hessian of the function.
//   */
//  using Hessian = Derivative<2>;
  ///@}

private:
  /**
   * Default constructor. It does nothing but it is needed for the serialization mechanism.
   */
  Domain() = default;

protected:
  Domain(std::shared_ptr<GridFuncType> func);

public:
  virtual ~Domain() = default;

  static std::shared_ptr<self_t>
  create(std::shared_ptr<GridFuncType> func)
  {
    return std::shared_ptr<self_t>(new self_t(func));
  }


  static std::shared_ptr<const self_t>
  const_create(std::shared_ptr<GridFuncType> func)
  {
    return create(func);
  }

  std::shared_ptr<GridFuncType> get_grid_function() const;

public:
  virtual std::unique_ptr<ElementHandler>
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


  void print_info(LogStream &out) const;

private:
  SharedPtrConstnessHandler<GridFunction<dim, space_dim>> grid_func_;

  friend class DomainElementBase<dim_, codim_, Domain<dim_, codim_>>;
  friend class DomainElementBase<dim_, codim_, const Domain<dim_, codim_>>;
  friend class DomainElement<dim_, codim_>;
  friend class ConstDomainElement<dim_, codim_>;

#ifdef SERIALIZATION
  /**
   * @name Functions needed for serialization
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar)
  {
    ar &grid_func_;
  }
  ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif

