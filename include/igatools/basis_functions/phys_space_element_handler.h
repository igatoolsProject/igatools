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

#ifndef PHYS_SPACE_ELEMENT_HANDLER_H_
#define PHYS_SPACE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/nurbs_element_handler.h>

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank, int codim,Transformation type_>
class PhysicalSpace;

/**
 * Element handler for an isogeometric space
 *
 * @ingroup serializable
 */
template<int dim_,int range_,int rank_,int codim_,Transformation type_>
class PhysSpaceElementHandler
  :
  public SpaceElementHandler<dim_,codim_,range_,rank_,type_>
{

  using PhysSpace = PhysicalSpace<dim_,range_,rank_,codim_,type_>;
  using RefSpace =  typename PhysSpace::RefSpace;
  using RefPhysSpaceElementHandler = typename PhysSpace::RefSpace::ElementHandler;
//    using PFCache = typename PhysSpace::PushForwardType;

  using ElementIterator = typename PhysSpace::ElementIterator;
  using ElementAccessor = typename PhysSpace::ElementAccessor;
  using PushFwd = typename PhysSpace::PushFwd;

  using base_t = SpaceElementHandler<dim_,codim_,range_,rank_,type_>;
  using self_t = PhysSpaceElementHandler<dim_,range_,rank_,codim_,type_>;

  using eval_pts_variant = SubElemVariants<Quadrature,dim_>;
  using topology_variant = TopologyVariants<dim_>;

public:
  static const int dim = dim_;

//    using PhysSpace::PushForwardType::type;

  /**
   * @name Constructors
   */
  ///@{
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  PhysSpaceElementHandler() = default;

protected:

  PhysSpaceElementHandler(std::shared_ptr<const PhysSpace> space);
  /**
   * Copy constructor. Not allowed to be used.
   */
  PhysSpaceElementHandler(const self_t &) = delete;

  /**
   * Move constructor. Not allowed to be used.
   */
  PhysSpaceElementHandler(self_t &&) = delete;

public:
  /**
   * Destructor.
   */
  virtual ~PhysSpaceElementHandler() = default;
  ///@}

  /**
   * Assignment operators.
   */
  ///@{
  /**
   * Copy assignment operator. Not allowed to be used.
   */
  self_t &operator=(const self_t &) = delete;

  /**
   * Move assignment operator. Not allowed to be used.
   */
  self_t &operator=(self_t &&) = delete;
  ///@}

  /**
   * @name Creators.
   */
  ///@{
  static std::shared_ptr<self_t> create(std::shared_ptr<const PhysSpace> space);
  ///@}

#if 0
  /**
   * Resets all the internal data in order to use the
   * quadrature scheme for the elements of the space with ID specified by
   * the input parameter <tt>elements_flat_id</tt>.
   */
  virtual void reset_selected_elements(
    const ValueFlags &flag,
    const eval_pts_variant &eval_points,
    const SafeSTLVector<int> &elements_flat_id) override final;


  virtual void init_cache(SpaceElement<dim_,codim_,range_,rank_,type_> &sp_elem,
                          const topology_variant &topology) override final;


  virtual void fill_cache(SpaceElement<dim_,codim_,range_,rank_,type_> &sp_elem,
                          const topology_variant &topology,
                          const int sub_elem_id) override final;
#endif

  void print_info(LogStream &out) const override final;

private:
  using RefElemHandler = SpaceElementHandler<RefSpace::dim,0,RefSpace::range,RefSpace::rank,Transformation::h_grad>;
  std::shared_ptr<RefElemHandler> ref_space_handler_;


//    using PushFwd = typename PhysSpace::PushForwardType;
//    PushFwd push_fwd_;

  using Map = typename PhysSpace::Map;
  Map mapping_;


  SafeSTLArray<ValueFlags, dim+1> flags_;

#if 0
  struct ResetDispatcher : boost::static_visitor<void>
  {
    ResetDispatcher(
      const ValueFlags flag_in,
      const SafeSTLVector<Index> &elements_flat_id,
      RefElemHandler &ref_space_handler,
      Map &mapping,
      SafeSTLArray<ValueFlags, dim+1> &flags)
      :
      flag_in_(flag_in),
      elements_flat_id_(elements_flat_id),
      ref_space_handler_(ref_space_handler),
      mapping_(mapping),
      flags_(flags)
    {};

    template<int sub_elem_dim>
    void operator()(const Quadrature<sub_elem_dim> &quad);

    const ValueFlags flag_in_;
    const SafeSTLVector<Index> &elements_flat_id_;
    RefElemHandler &ref_space_handler_;
    Map &mapping_;
    SafeSTLArray<ValueFlags, dim+1> &flags_;
  };


  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(
      const SafeSTLArray<ValueFlags, dim+1> &flags,
      RefElemHandler &ref_space_handler,
      Map &mapping,
      PhysicalSpaceElement<dim_,range_,rank_,codim_,type_> &phys_elem)
      :
      flags_(flags),
      ref_space_handler_(ref_space_handler),
      mapping_(mapping),
      phys_elem_(phys_elem)
    {};

    template<int sub_elem_dim>
    void operator()(const Topology<sub_elem_dim> &topology);

    const SafeSTLArray<ValueFlags, dim+1> &flags_;
    RefElemHandler &ref_space_handler_;
    Map &mapping_;
    PhysicalSpaceElement<dim_,range_,rank_,codim_,type_> &phys_elem_;
  };



  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(
      const int sub_elem_id,
      RefElemHandler &ref_space_handler,
      Map &mapping,
      PhysicalSpaceElement<dim_,range_,rank_,codim_,type_> &phys_elem)
      :
      sub_elem_id_(sub_elem_id),
      ref_space_handler_(ref_space_handler),
      mapping_(mapping),
      phys_elem_(phys_elem)
    {};

    template<int sub_elem_dim>
    void operator()(const Topology<sub_elem_dim> &topology);

    const int sub_elem_id_;
    RefElemHandler &ref_space_handler_;
    Map &mapping_;
    PhysicalSpaceElement<dim_,range_,rank_,codim_,type_> &phys_elem_;
  };

#endif


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
    ar &boost::serialization::make_nvp("PhysSpaceElementHandler_base_t",
                                       boost::serialization::base_object<base_t>(*this));

    ar.template register_type<BSplineElementHandler<dim_,range_,rank_> >();
    ar.template register_type<NURBSElementHandler<dim_,range_,rank_> >();
    ar &boost::serialization::make_nvp("ref_space_handler_",ref_space_handler_);


//        ar &boost::serialization::make_nvp("push_fwd_",push_fwd_);


    ar &boost::serialization::make_nvp("flags_",flags_);
  }
  ///@}
#endif

};


IGA_NAMESPACE_CLOSE

#endif
