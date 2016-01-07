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
#include <igatools/geometry/grid.h>
#include <igatools/geometry/grid_element.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Grid caches manager
 *
 * @ingroup handlers
 */
template <int dim>
class GridHandler
{
private:
  using self_t = GridHandler<dim>;

public:
  using GridType = const Grid<dim>;


protected:
  using ElementIterator = typename GridType::ElementIterator;
  using ElementAccessor = typename GridType::ElementAccessor;

public:
  using Flags = grid_element::Flags;
protected:
  using FlagsArray = SafeSTLArray<Flags, dim+1>;
  using topology_variant = TopologyVariants<dim>;
  template<int k>
  using ConstQuad = const Quadrature<k>;
  using eval_pts_variant = SubElemPtrVariants<ConstQuad,dim>;


  /**
   * @name Constructors
   */
  ///@{
protected:

  GridHandler() = delete;

public:
  /**
   * Constructor.
   */
  GridHandler(const std::shared_ptr<GridType> &grid);

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

  void set_element_flags(const Flags &flag);


  template <int sdim>
  void init_cache(ElementAccessor &elem,
                  const std::shared_ptr<const Quadrature<sdim>> &quad) const;

  template <int sdim>
  void init_cache(ElementIterator &elem,
                  const std::shared_ptr<const Quadrature<sdim>> &quad) const;

  void init_cache(ElementAccessor &elem,
                  const eval_pts_variant &quad) const;



  void init_element_cache(ElementIterator &elem,
                          const std::shared_ptr<const Quadrature<dim>> &quad) const;

  void init_face_cache(ElementIterator &elem,
                       const std::shared_ptr<const Quadrature<(dim > 0) ? dim-1 : 0>> &quad) const;

  template <int sdim>
  void fill_cache(ElementAccessor &elem, const int s_id) const;

  template <int sdim>
  void fill_cache(ElementIterator &elem, const int s_id) const;

  void fill_cache(const topology_variant &sdim,
                  ElementAccessor &elem,
                  const int s_id) const;





  void fill_element_cache(ElementIterator &elem) const;


  void fill_face_cache(ElementIterator &elem, const int s_id) const;
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


  struct SetFlagsDispatcher : boost::static_visitor<void>
  {
    SetFlagsDispatcher(const Flags flag, self_t &grid_handler);

    template<int sdim>
    void operator()(const Topology<sdim> &s_el);

    const Flags flag_;
    self_t &grid_handler_;
  };


  struct InitCacheDispatcher : boost::static_visitor<void>
  {
    InitCacheDispatcher(const self_t &grid_handler,
                        ElementAccessor &elem);

    template<int sdim>
    void operator()(const std::shared_ptr<const Quadrature<sdim>> &quad);

    const self_t &grid_handler_;
    ElementAccessor &elem_;

  };


  struct FillCacheDispatcher : boost::static_visitor<void>
  {
    FillCacheDispatcher(const self_t &grid_handler,
                        ElementAccessor &elem,
                        const int s_id);

    template<int sdim>
    void operator()(const Topology<sdim> &);

    const self_t &grid_handler_;
    ElementAccessor &elem_;
    int s_id_;
  };

};

IGA_NAMESPACE_CLOSE

#endif /* GRID_HANDLER_H_ */
