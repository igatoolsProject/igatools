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
//#include <igatools/geometry/domain_handler.h>

//#include <igatools/base/tensor.h>
//#include <igatools/geometry/unit_element.h>
//#include <igatools/utils/value_vector.h>
//#include <igatools/base/quadrature.h>
//#include <igatools/geometry/grid.h>
//#include <igatools/geometry/grid_iterator.h>

IGA_NAMESPACE_OPEN

template <int, int, int, int> class FunctionElement;
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

#if 0
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
#endif

  std::shared_ptr<const DomainType> get_domain() const
  {
    return domain_.get_ptr_const_data();
  }



public:
  virtual std::unique_ptr<ElementHandler>
  create_cache_handler() const;

#if 0
  std::unique_ptr<ElementAccessor>
  create_element(const ListIt &index, const PropId &prop) const;
#endif

  std::unique_ptr<ElementAccessor>
  create_element_begin(const PropId &prop) const;

  std::unique_ptr<ElementAccessor>
  create_element_end(const PropId &prop) const;

public:
  ///@name Iterating of grid elements
  ///@{
#if 0
  /**
   * This function returns a element iterator to the first element of the patch.
   */
  ElementIterator begin(const PropId &property = ElementProperties::active);

  /**
   * This function returns a element iterator to one-pass the end of patch.
   */
  ElementIterator end(const PropId &property = ElementProperties::active);
#endif

  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  ElementIterator begin(const PropId &property = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  ElementIterator end(const PropId &property = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to the first element of the patch.
   */
  ElementIterator cbegin(const PropId &property = ElementProperties::active) const;

  /**
   * This function returns a element (const) iterator to one-pass the end of patch.
   */
  ElementIterator cend(const PropId &property = ElementProperties::active) const;
  ///@}

  virtual void print_info(LogStream &out) const = 0;


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




  template <int sdim>
  using SubFunc = Function<(sdim>0)?sdim:0,codim+(dim-((sdim>0)?sdim:0)),range,rank>;

  virtual
  std::shared_ptr<const SubFunc<(dim>0)?dim-1:0> >
  get_sub_function(const int s_id,const std::shared_ptr<const Grid<(dim>0)?dim-1:0>> &sub_grid = nullptr) const
  {
    Assert(dim >= 1,ExcMessage("The topological dimension must be >=1"));
#if 0
    const int sdim = (dim>0)?dim-1:0;

    using SubGridElemMap = typename Grid<dim>::template SubGridMap<sdim>;
    SubGridElemMap sub_grid_elem_map;
    using SubGrid = Grid<sdim>;

    std::shared_ptr<const SubGrid> s_grid;

    if (sub_grid == nullptr)
    {
      const auto grid = domain_->get_grid_function()->get_grid();
      s_grid = grid->template get_sub_grid<sdim>(s_id,sub_grid_elem_map);
    }
    else
    {
      s_grid = sub_grid;
    }

    auto bndry_domain = domain_->get_sub_domain(s_id,sub_grid_elem_map,s_grid);
#endif


    AssertThrow(false,ExcMessage("This function must be implemented by a derived class."));
    AssertThrow(false,ExcNotImplemented());
    return nullptr;
  }

protected:
  SharedPtrConstnessHandler<DomainType> domain_;

  /**
   * Name associated to the object instance.
   */
  std::string name_;

  /**
   * Unique identifier associated to each object instance.
   */
  Index object_id_;


private:
  friend class FunctionElement<dim_, codim_, range_, rank_>;




#ifdef MESH_REFINEMENT
protected:
  std::shared_ptr<const self_t> function_previous_refinement_;

public:
  const std::shared_ptr<const self_t> &get_function_previous_refinement() const;

private:
  /**
   * Rebuild the internal state of the object after an insert_knots() function is invoked.
   *
   * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
   * @note This function is connected to the Grid's signal for the refinement, and
   * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
   *
   * @ingroup h_refinement
   */
  virtual void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid)
  {
    AssertThrow(false,ExcMessage("This function must be implemented in a derived class."));
  }

public:

  /**
   *  Connect a slot (i.e. a function pointer) to the refinement signals
   *  which will be
   *  emitted whenever a insert_knots() function is called by the underlying
   *  a Grid member.
   */
  boost::signals2::connection
  connect_insert_knots(const typename Grid<dim_>::SignalInsertKnotsSlot &subscriber);

  void create_connection_for_insert_knots(const std::shared_ptr<self_t> &function);

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
  serialize(Archive &ar);
  ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif
