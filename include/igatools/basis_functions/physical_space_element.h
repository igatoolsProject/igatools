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


#ifndef PHYSICAL_SPACE_ELEMENT_H
#define PHYSICAL_SPACE_ELEMENT_H

#include <igatools/base/config.h>

#include <igatools/base/quadrature.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/nurbs_element.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/geometry/push_forward.h>



IGA_NAMESPACE_OPEN

template <class Accessor> class GridIterator;

/**
 *
 * @ingroup elements
 * @ingroup serializable
 */
template<int dim_,int range_,int rank_,int codim_,Transformation type_ = Transformation::h_grad>
class PhysicalSpaceElement
  :
  public SpaceElement<dim_,codim_,range_,rank_,type_>
{
  using self_t = PhysicalSpaceElement<dim_,range_,rank_,codim_,type_>;
  using parent_t = SpaceElement<dim_,codim_,range_,rank_,type_>;

public :
  using PhysSpace = PhysicalSpace<dim_,range_,rank_,codim_,type_>;
  /** Type required by the GridIterator templated iterator */
  using ContainerType = const PhysSpace;

  using Space = PhysSpace;
  using RefSpace = typename PhysSpace::RefSpace;
  using PushFwd = typename PhysSpace::PushFwd;
//    using RefElemAccessor = SpaceElement<RefSpace::dim,0,RefSpace::range,RefSpace::rank,Transformation::h_grad>;
  using RefElemAccessor = ReferenceElement<RefSpace::dim,RefSpace::range,RefSpace::rank>;

  using PhysDomain = Domain<dim_, codim_>;
  using PhysDomainElem = DomainElement<dim_, codim_>;

  static const auto dim = PushFwd::dim;
  static const auto space_dim = PushFwd::space_dim;
  static const auto codim = PushFwd::codim;
  static const auto type = PushFwd::type;

  using PhysPoint = typename Space::Point;

  using Grid = Grid<dim>;
  using IndexType = typename Grid::IndexType;
  using List = typename Grid::List;
  using ListIt = typename Grid::ListIt;


  /**
   * @name Constructors
   */
  ///@{
public:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  PhysicalSpaceElement() = default;

  PhysicalSpaceElement(const std::shared_ptr<ContainerType> space,
                       const ListIt &index,
                       const PropId &prop = ElementProperties::active);


  /**
   * Copy constructor. Not allowed to be used.
   */
  PhysicalSpaceElement(const self_t &in) = delete;


  /**
   * Move constructor.
   */
  PhysicalSpaceElement(self_t &&in) = default;

  /**
   * Destructor.
   */
  virtual ~PhysicalSpaceElement() = default;

  ///@}

  /**
   * @name Assignment operators
   */
  ///@{
  /**
   * Copy assignment operator. Not allowed to be used.
   */
  self_t &
  operator=(const self_t &in) = delete;

  /**
   * Move assignment operator.
   */
  self_t &
  operator=(self_t &&in) = default;

  ///@}



  /**
   * @name Getting quantities that are geometry-related
   */
  ///@{
  /**
   * Returns the <tt>k</tt> dimensional j-th sub-element measure
   * multiplied by the weights of the quadrature.
   */
  template <int k>
  const ValueVector<Real> &get_w_measures(const int j) const
  {
    return phys_domain_element_->template get_w_measures<k>(j);
  }

  /**
   * Returns the gradient determinant of the map at the dilated quadrature points.
   */
  template <int k>
  const ValueVector<Real> &get_measures(const int j) const
  {
    return phys_domain_element_->template get_measures<k>(j);
  }

  const ValueVector<Real> &get_element_w_measures() const;

  template <int k = dim>
  const ValueVector<PhysPoint> &get_points(const int j = 0) const;

  const ValueVector<PhysPoint> &get_element_points() const;

  template<int sub_dim>
  const ValueVector<Points<space_dim> > &
  get_boundary_normals(const int s_id) const
  {
    return phys_domain_element_->template get_boundary_normals<sub_dim>(s_id);
  }


  /**
   * Prints internal information about the BSplineElementAccessor.
   * Its main use is for testing and debugging.
   */
  void print_info(LogStream &out) const override final;

  void print_cache_info(LogStream &out) const;

#if 0
  /**
   * @name Functions for the basis evaluation without the use of the cache.
   */
  ///@{
  /**
   * Returns a ValueTable with the quantity specified by the template parameter @p ValueType,
   * computed for all local basis function,
   * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
   * @note This function does not use the cache and therefore can be called any time without
   * needing to pre-call init_cache()/fill_cache().
   * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
   * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
   */
  template <class ValueType>
  decltype(auto)
  evaluate_basis_at_points(
    const Quadrature<dim_> &points,
    const std::string &dofs_property)
  {
    auto elem_handler = typename Space::ElementHandler::create(this->space_);

    elem_handler->reset_one_element(ValueType::flag,points,this->get_flat_index());
    elem_handler->template init_cache<dim>(*this);
    elem_handler->template fill_cache<dim>(*this,0);

    return this->template get_basis<ValueType,dim>(0,dofs_property);
  }
#endif
  ///@}

public:

  /**
   * Return a const reference of the reference space element.
   */
  const RefElemAccessor &get_ref_space_element() const;

  /**
   * Return a non-const reference of the reference space element.
   */
  RefElemAccessor &get_ref_space_element();


  /**
   * Return a const reference of the DomainElement.
   */
  const PhysDomainElem &get_physical_domain_element() const;

  /**
   * Return a non-const reference of the DomainElement.
   */
  PhysDomainElem &get_physical_domain_element();

public:
  using parent_t::get_num_basis;



  /** Returns the index of the element. */
  IndexType get_index() const;


  /** Return the cartesian grid from which the element belongs.*/
  const std::shared_ptr<const Grid<dim>> get_grid() const;


  virtual typename List::iterator &operator++() override final
  {
    ++(*phys_domain_element_);
    return ++(*ref_space_element_);
    Assert(false,ExcNotImplemented());
  }


#if 0
  /**
   * For a given flags input argument identifies the face quantities and
   * returns a new ValueFlags variable containing only face quantities.
   * The output flags does not contain the word face.
   */
  ValueFlags get_face_flags(const ValueFlags fill_flag) const ;

#endif

#if 0
  /** @name Functions/operators for moving the element in the Grid.*/
  ///@{
  /**
   * Sets the index of the element using the flatten representation.
   * @note This function also updates the index for the tensor representation.
   * @warning This may be a dangerous function, be careful when using it
   * as it is easy to use incorrectly. Only use it if you know what you
   * are doing.
   */
  void move_to(const Index flat_index) override final;
  ///@}
#endif

protected:

  /*
      bool operator==(const PhysicalSpaceElement<PhysSpace> &a) const;
      bool operator!=(const PhysicalSpaceElement<PhysSpace> &a) const;
      bool operator<(const PhysicalSpaceElement<PhysSpace> &a) const;
      bool operator>(const PhysicalSpaceElement<PhysSpace> &a) const;
  //*/

#if 0
  /**
   * This function returns the ValueFlags needed to be passed to the ReferenceSpacePhysicalAccessor
   * in order to compute the quantities specified by the input argument
   * @p fill_flag (i.e. the ValueFlags that refers to the PhysicalSpaceElement).
   */
  ValueFlags get_reference_space_accessor_fill_flags(const ValueFlags fill_flag) const;

  /**
   * This function returns the ValueFlags needed to be passed to the PushForwardAccessor
   * in order to compute the quantities specified by the input argument
   * @p fill_flag (i.e. the ValueFlags that refers to the PhysicalSpaceElement).
   */
  ValueFlags get_push_forward_accessor_fill_flags(const ValueFlags fill_flag) const;

#endif
  /**
   * Performs a copy of the input @p element.
   * The type of copy (deep or shallow) is specified by the input parameter @p copy_policy.
   */
  void copy_from(const self_t &element,
                 const CopyPolicy &copy_policy);


private:
  template <class Accessor> friend class GridIteratorBase;
  template <int,int,int,int,Transformation> friend class PhysSpaceElementHandler;

  std::unique_ptr<typename RefElemAccessor::parent_t> ref_space_element_;

  std::unique_ptr<PhysDomainElem> phys_domain_element_;


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
    ar &boost::serialization::make_nvp("PhysicalSpaceElement_base_t",
                                       boost::serialization::base_object<parent_t>(*this));

    ar.template register_type<BSplineElement<dim_,range_,rank_> >();
    ar.template register_type<NURBSElement<dim_,range_,rank_> >();
    ar &boost::serialization::make_nvp("ref_space_element_",ref_space_element_);


    ar &boost::serialization::make_nvp("phys_domain_element_",phys_domain_element_);
  }
  ///@}
#endif

};


IGA_NAMESPACE_CLOSE

#endif
