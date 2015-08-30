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

#ifndef __GRID_CACHE_HANDLER_H_
#define __GRID_CACHE_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/values_cache.h>
#include <igatools/base/tuple_utils.h>
#include <igatools/base/quadrature.h>
#include <igatools/utils/tensor_product_array.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/grid_element.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Grid caches manager
 *
 * @ingroup serializable
 */
template <int dim>
class GridElementHandler
{
private:
  using self_t = GridElementHandler<dim>;

public:
  using GridType = const CartesianGrid<dim>;


protected:
  // using ElementIterator = typename GridType::ElementIterator;
  // using ElementAccessor = typename GridType::ElementAccessor;
  using ElementIterator = typename GridType::ElementConstIterator;
  using ElementAccessor = typename GridType::ConstElementAccessor;

public:
  using Flags = typename ElementAccessor::Flags;
protected:
  using FlagsArray = SafeSTLArray<Flags, dim+1>;
  using topology_variant = TopologyVariants<dim>;

public:
  /**
   * @name Creators.
   */
  ///@{
  static std::shared_ptr<self_t> create(std::shared_ptr<GridType> grid);
  ///@}

  /**
   * @name Constructors
   */
  ///@{
protected:
  /**
   * Default constructor. It does nothing but it is needed for the
   * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   * mechanism of the Function class.
   */
  GridElementHandler() = default;

public:
  /**
   * Constructor.
   */
  GridElementHandler(std::shared_ptr<GridType> grid);

  /**
   * Copy constructor.
   */
  GridElementHandler(const self_t &) = default;

  /**
   * Move constructor.
   */
  GridElementHandler(self_t &&) = default;

  /**
   * Destructor.
   */
  ~GridElementHandler() = default;
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

public:
  /**
   * @name Functions for the cache's reset/init/fill mechanism.
   */
  ///@{
  template<int sdim>
  void set_flags(const Flags &flag);

  void set_flags(const topology_variant &sdim,
		  const Flags &flag);

  template <int sdim>
  void init_cache(ElementAccessor &elem,
                  std::shared_ptr<const Quadrature<sdim>> quad) const;

  template <int sdim>
  void init_cache(ElementIterator &elem,
                  std::shared_ptr<const Quadrature<sdim>> quad) const
  {
    init_cache<sdim>(*elem, quad);
  }


  void init_element_cache(ElementIterator &elem,
                          std::shared_ptr<const Quadrature<dim>> quad) const
  {
    init_cache<dim>(*elem, quad);
  }

  void init_face_cache(ElementIterator &elem,
                       std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> quad) const
  {
    Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
    init_cache<(dim > 0) ? dim-1 : 0>(*elem, quad);
  }

  template <int sdim>
  void fill_cache(ElementAccessor &elem, const int s_id) const;


  template <int sdim>
  void fill_cache(ElementIterator &elem, const int s_id) const
  {
    fill_cache<sdim>(*elem, s_id);
  }



  void fill_element_cache(ElementIterator &elem) const
  {
    fill_cache<dim>(*elem,0);
  }


  void fill_face_cache(ElementIterator &elem, const int s_id) const
  {
    Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
    fill_cache<(dim > 0) ? dim-1 : 0>(*elem,s_id);
  }
  ///@}


public:

  /**
   * Function for printing some internal information.
   * Its use is mostly intended for debugging and testing purposes.
   */
  void print_info(LogStream &out) const;


  /**
   * Returns the grid upon which the object is built.
   */
  std::shared_ptr<const GridType> get_grid() const;

private:
  std::shared_ptr<GridType> grid_;

  FlagsArray flags_;


private:
  struct SetFlagsDispatcher : boost::static_visitor<void>
    {
  	  SetFlagsDispatcher(const Flags flag, FlagsArray &flags)
          		:
          			flag_(flag),
  					flags_(flags)
  					{}

  	  template<int sdim>
  	  void operator()(const Topology<sdim> &)
  	  {
  		 // grid_handler_.set_flags<sdim>(flag_)
  		  flags_[sdim] = flag_;
  	  }

  	  const Flags flag_;
  	  FlagsArray &flags_;
    };

private:

#ifdef SERIALIZATION
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class boost::serialization::access;

  template<class Archive>
  void
  serialize(Archive &ar, const unsigned int version);
  ///@}
#endif // SERIALIZATION

};

IGA_NAMESPACE_CLOSE

#endif /* GRID_ELEMENT_HANDLER_H_ */
