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

#ifndef __GRID_HANDLER_H_
#define __GRID_HANDLER_H_

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
class GridHandler
{
private:
  using self_t = GridHandler<dim>;

public:
  using GridType = const Grid<dim>;


protected:
  // using ElementIterator = typename GridType::ElementIterator;
  // using ElementAccessor = typename GridType::ElementAccessor;
  using ElementConstIterator = typename GridType::ElementConstIterator;
  using ConstElementAccessor = typename GridType::ConstElementAccessor;

public:
  using Flags = typename ConstElementAccessor::Flags;
protected:
  using FlagsArray = SafeSTLArray<Flags, dim+1>;
  using topology_variant = TopologyVariants<dim>;
  template<int k>
  using ConstQuad = const Quadrature<k>;
  using eval_pts_variant = SubElemPtrVariants<ConstQuad,dim>;

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
  GridHandler() = default;

public:
  /**
   * Constructor.
   */
  GridHandler(std::shared_ptr<GridType> grid);

  /**
   * Copy constructor.
   */
  GridHandler(const self_t &) = default;

  /**
   * Move constructor.
   */
  GridHandler(self_t &&) = default;

  /**
   * Destructor.
   */
  ~GridHandler() = default;
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
  void init_cache(ConstElementAccessor &elem,
                  std::shared_ptr<const Quadrature<sdim>> quad) const;

  void init_cache(ConstElementAccessor &elem,
                  const eval_pts_variant &quad) const;

  template <int sdim>
  void init_cache(ElementConstIterator &elem,
                  std::shared_ptr<const Quadrature<sdim>> quad) const
  {
    init_cache<sdim>(*elem, quad);
  }


  void init_element_cache(ElementConstIterator &elem,
                          std::shared_ptr<const Quadrature<dim>> quad) const
  {
    init_cache<dim>(*elem, quad);
  }

  void init_face_cache(ElementConstIterator &elem,
                       std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> quad) const
  {
    Assert(dim > 0,ExcMessage("No face defined for element with topological dimension 0."));
    init_cache<(dim > 0) ? dim-1 : 0>(*elem, quad);
  }

  template <int sdim>
  void fill_cache(ConstElementAccessor &elem, const int s_id) const;

  void fill_cache(const topology_variant &sdim,
                  ConstElementAccessor &elem,
                  const int s_id) const;


  template <int sdim>
  void fill_cache(ElementConstIterator &elem, const int s_id) const
  {
    fill_cache<sdim>(*elem, s_id);
  }



  void fill_element_cache(ElementConstIterator &elem) const
  {
    fill_cache<dim>(*elem,0);
  }


  void fill_face_cache(ElementConstIterator &elem, const int s_id) const
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
  std::shared_ptr<GridType> get_grid() const;

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


  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(self_t const *grid_handler,
                        ConstElementAccessor &elem)
      :
      grid_handler_(grid_handler),
      elem_(elem)
    {}

    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad)
    {
      grid_handler_->template init_cache<sdim>(elem_, quad);
    }

    self_t const *grid_handler_;
    ConstElementAccessor &elem_;

  };




  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(self_t const *grid_handler,
                        ConstElementAccessor &elem,
                        const int s_id)
      :
      grid_handler_(grid_handler),
      elem_(elem),
      s_id_(s_id)
    {}

    template<int sdim>
    void operator()(const Topology<sdim> &)
    {
      grid_handler_->template fill_cache<sdim>(elem_, s_id_);
    }

    self_t const *grid_handler_;
    ConstElementAccessor &elem_;
    int s_id_;

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
