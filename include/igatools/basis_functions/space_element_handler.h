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


#ifndef SPACE_ELEMENT_HANDLER_H_
#define SPACE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/space.h>
#include <igatools/base/tuple_utils.h>
#include <igatools/basis_functions/space_element.h>

IGA_NAMESPACE_OPEN



/**
 *
 */
template <int dim,int codim,int range,int rank,Transformation type>
class SpaceElementHandler
{
private:
  using self_t = SpaceElementHandler<dim,codim,range,rank,type>;

public:
  using Sp = Space<dim,codim,range,rank,type>;
  using ElementAccessor = typename Sp::ElementAccessor;
  using ElementIterator = typename Sp::ElementIterator;


  using IndexType = typename Sp::IndexType;


private:
  using eval_pts_variant = QuadVariants<dim>;
  using topology_variant = TopologyVariants<dim>;


protected:
  /** @name Constructors */
  ///@{
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism.
   */
  SpaceElementHandler() = default;


  SpaceElementHandler(std::shared_ptr<const Sp> space);

  /**
   * Copy constructor. Not allowed to be used.
   */
  SpaceElementHandler(const self_t &elem_handler) = delete;

  /**
   * Move constructor. Not allowed to be used.
   */
  SpaceElementHandler(self_t &&elem_handler) = delete;

public:

  /**
   * Destructor.
   */
  virtual ~SpaceElementHandler() = default;

  ///@}


  virtual void set_flags_impl(const topology_variant &topology,
                              const typename space_element::Flags &flag) = 0;

  template<int sdim>
  void set_flags(const typename space_element::Flags &flag)
  {
    this->set_flags_impl(Topology<sdim>(),flag);
  }


  template <int sdim>
  void init_cache(ElementAccessor &elem,
                  const std::shared_ptr<const Quadrature<sdim>> &quad) const
  {
    this->init_cache_impl(elem,quad);
  }

  template <int sdim>
  void init_cache(ElementIterator &elem,
                  std::shared_ptr<const Quadrature<sdim>> quad) const
  {
    init_cache<sdim>(*elem, quad);
  }

  void init_element_cache(ElementAccessor &elem,
                          std::shared_ptr<const Quadrature<dim>> quad) const
  {
    init_cache<dim>(elem, quad);
  }

  void init_element_cache(ElementIterator &elem,
                          std::shared_ptr<const Quadrature<dim>> quad) const
  {
    init_element_cache(*elem, quad);
  }

  void init_face_cache(ElementAccessor &elem,
                       std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> quad) const
  {
    Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
    init_cache<(dim > 0) ? dim-1 : 0>(elem, quad);
  }

  void init_face_cache(ElementIterator &elem,
                       std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> quad) const
  {
    init_face_cache(*elem, quad);
  }


  template <int sdim>
  void fill_cache(ElementAccessor &elem, const int s_id) const
  {
    Assert(s_id >= 0 && s_id < UnitElement<dim>::template num_elem<sdim>(),
           ExcIndexRange(s_id,0,UnitElement<dim>::template num_elem<sdim>()));
    this->fill_cache_impl(elem,Topology<sdim>(),s_id);
  }


  template <int sdim>
  void fill_cache(ElementIterator &elem, const int s_id) const
  {
    fill_cache<sdim>(*elem, s_id);
  }

  void fill_element_cache(ElementAccessor &elem) const
  {
    fill_cache<dim>(elem,0);
  }

  void fill_element_cache(ElementIterator &elem) const
  {
    fill_element_cache(*elem);
  }


  void fill_face_cache(ElementAccessor &elem, const int s_id) const
  {
    Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
    fill_cache<(dim > 0) ? dim-1 : 0>(elem,s_id);
  }

  void fill_face_cache(ElementIterator &elem, const int s_id) const
  {
    fill_face_cache(*elem,s_id);
  }


public:


#if 0
  /**
   * Resets all the internal data in order to use the
   * same quadrature scheme for each active element of the space.
   */
  void reset(const ValueFlags &flag, const eval_pts_variant &quad);


  /**
   * Resets all the internal data in order to use the
   * same quadrature scheme for the elements of the space with ID specified by
   * the input parameter <tt>elements_flat_id</tt>.
   *
   * @note This function is pure virtual and must be implemented in the class that are derived
   * from SpaceElementHandler.
   */
  virtual void reset_selected_elements(
    const ValueFlags &flag,
    const eval_pts_variant &eval_points,
    const SafeSTLVector<IndexType> &elements_id) = 0;


  virtual void init_cache(SpaceElement<dim,codim,range,rank,type> &elem,
                          const topology_variant &topology) = 0;

  template <int sub_elem_dim>
  void init_cache(SpaceElement<dim,codim,range,rank,type> &elem);

  /**
   * Allocates the space in the cache of ElementAccessor <tt>element</tt>
   * necessary for the given quadrature and flag combination.
   * It also fills the invariant (not changing) members of
   * the cache.
   */
  void init_element_cache(ElementAccessor &elem);

  void init_element_cache(ElementIterator &elem);

  void init_face_cache(ElementAccessor &elem);

  void init_face_cache(ElementIterator &elem);


  virtual void fill_cache(
    SpaceElement<dim,codim,range,rank,type> &elem,
    const topology_variant &topology,
    const int sub_elem_id) = 0;

  template<int sub_elem_dim>
  void fill_cache(ElementAccessor &elem, const int sub_elem_id);

  void fill_element_cache(ElementAccessor &elem);

  void fill_element_cache(ElementIterator &elem);

  void fill_face_cache(ElementAccessor &elem, const int face_id);

  void fill_face_cache(ElementIterator &elem, const int face_id);
#endif

  virtual void print_info(LogStream &out) const = 0;


  std::shared_ptr<const Sp> get_space() const;

private:
  std::shared_ptr<const Sp> space_;



  virtual void init_cache_impl(ElementAccessor &elem,
                               const eval_pts_variant &quad) const = 0;

  virtual void fill_cache_impl(ElementAccessor &elem,
                               const topology_variant &topology,
                               const int s_id) const = 0;


protected:
  SafeSTLArray<typename space_element::Flags, dim + 1> flags_;

};

IGA_NAMESPACE_CLOSE

#endif // SPACE_ELEMENT_HANDLER_H_
