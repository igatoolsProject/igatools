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

#ifndef __FUNCTION_H_
#define __FUNCTION_H_

#include <igatools/base/config.h>
#include <igatools/geometry/domain.h>
#include <igatools/geometry/domain_handler.h>

//#include <igatools/base/tensor.h>
//#include <igatools/geometry/unit_element.h>
//#include <igatools/utils/value_vector.h>
//#include <igatools/base/quadrature.h>
//#include <igatools/geometry/grid.h>
//#include <igatools/geometry/grid_iterator.h>

IGA_NAMESPACE_OPEN

template <int, int, int, int, class> class FunctionElementBase;
template <int, int, int, int> class FunctionElement;
template <int, int, int, int> class ConstFunctionElement;
template <int, int, int, int> class FunctionHandler;

/**
 * Function Class
 *
 * @ingroup serializable
 */
template<int dim_, int codim_ = 0, int range_ = 1, int rank_ = 1>
class Function :
  public std::enable_shared_from_this<Function<dim_,codim_,range_,rank_> >
{
private:
  using self_t = Function<dim_, codim_, range_, rank_>;

public:
  static const int space_dim = dim_ + codim_;
  static const int dim       = dim_;
  static const int codim     = codim_;
  static const int range     = range_;
  static const int rank      = rank_;

  using DomainType = Domain<dim_, codim_>;

  using ElementAccessor = FunctionElement<dim_, codim_, range_, rank_>;
  using ElementIterator = GridIterator<ElementAccessor>;
  using ConstElementAccessor = ConstFunctionElement<dim_, codim_, range_, rank_>;
  using ElementConstIterator = GridIterator<ConstElementAccessor>;

  using ElementHandler = FunctionHandler<dim_, codim_, range_, rank_>;

  using List = typename DomainType::List;
  using ListIt = typename DomainType::ListIt;

  /** Types for the input/output evaluation arguments */
  ///@{
  /**
   * Type for the input argument of the function.
   */
  using Point = Points<space_dim>;

  /**
   * Type for the return of the function.
   */
  using Value = Values<space_dim, range_, rank_>;

  /**
   * Type for the derivative of the function.
   */
  template <int order>
  using Derivative = Derivatives<space_dim, range_, rank_, order>;

  /**
   * Type for the gradient of the function.
   */
  using Gradient = Derivative<1>;

  /**
   * Type for the hessian of the function.
   */
  using Hessian = Derivative<2>;

  /**
   * Type for the divergence of function.
   */
  using Div = Values<space_dim, space_dim, rank_-1>;
  ///@}

  /** @name Constructors and destructor. */
  ///@{

public:
  /**
   * Default constructor. It does nothing but it is needed for the serialization
   * mechanism.
   */
  Function() = default;

protected:
  /** Constructor */
  Function(const SharedPtrConstnessHandler<DomainType> &domain);


public:
  /** Destructor */
  virtual ~Function() = default;
  ///@}

  static std::shared_ptr<self_t>
  create(std::shared_ptr<DomainType> domain)
  {
    return std::shared_ptr<self_t>(new
                                   self_t(SharedPtrConstnessHandler<DomainType>(domain)));
  }


  static std::shared_ptr<const self_t>
  const_create(std::shared_ptr<const DomainType> domain)
  {
    return std::shared_ptr<self_t>(new self_t(
                                     SharedPtrConstnessHandler<DomainType>(domain)));
  }


  std::shared_ptr<const DomainType> get_domain() const
  {
    return domain_.get_ptr_const_data();
  }



public:
  virtual std::unique_ptr<ElementHandler>
  create_cache_handler() const;

  std::unique_ptr<ConstElementAccessor>
  create_element(const ListIt &index, const PropId &prop) const;

  std::unique_ptr<ElementAccessor>
  create_element(const ListIt &index, const PropId &prop);

public:
  ///@name Iterating of grid elements
  ///@{
  /**
   * This function returns a element iterator to the first element of the patch.
   */
  ElementIterator begin(const PropId &property = ElementProperties::active);

  /**
   * This function returns a element iterator to one-pass the end of patch.
   */
  ElementIterator end(const PropId &property = ElementProperties::active);

  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  ElementConstIterator begin(const PropId &property = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  ElementConstIterator end(const PropId &property = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  ElementConstIterator cbegin(const PropId &property = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  ElementConstIterator cend(const PropId &property = ElementProperties::active) const;
  ///@}

  virtual void print_info(LogStream &out) const
  {
    Assert(false, ExcNotImplemented());
  }


  /**
   * Get the name associated to the object instance.
   */
  const std::string &get_name() const;

  /**
   * Set the name associated to the object instance.
   */
  void set_name(const std::string &name);

  /**
   * Returns the unique identifier associated to each object instance.
   */
  Index get_object_id() const;

private:
  SharedPtrConstnessHandler<DomainType> domain_;

  /**
   * Name associated to the object instance.
   */
  std::string name_;

  /**
   * Unique identifier associated to each object instance.
   */
  Index object_id_;


  friend class FunctionElementBase<dim_, codim_, range_, rank_, Function<dim_, codim_, range_, rank_>>;
  friend class FunctionElementBase<dim_, codim_, range_, rank_, const Function<dim_, codim_, range_, rank_>>;
  friend class FunctionElement<dim_, codim_, range_, rank_>;
  friend class ConstFunctionElement<dim_, codim_, range_, rank_>;



#ifdef MESH_REFINEMENT
protected:
  std::shared_ptr<const self_t> function_previous_refinement_;

public:
  const std::shared_ptr<const self_t> &get_function_previous_refinement() const
  {
    return function_previous_refinement_;
  }
#endif // MESH_REFINEMENT

#ifdef SERIALIZATION
private:
  /**
   * @name Functions needed for the serialization
   */
  ///@{
  friend class cereal::access;

  template<class Archive>
  void
  serialize(Archive &ar)
  {
    ar &make_nvp("domain_",domain_);
    ar &make_nvp("name_",name_);
    ar &make_nvp("object_id_",object_id_);

#ifdef MESH_REFINEMENT
    auto tmp = std::const_pointer_cast<self_t>(function_previous_refinement_);
    ar &make_nvp("function_previous_refinement_",tmp);
    Assert(tmp != nullptr,ExcNullPtr());

    function_previous_refinement_ = tmp;
#endif // MESH_REFINEMENT
  }
  ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif
