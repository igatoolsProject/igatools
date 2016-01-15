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

#ifndef __SUB_GRID_FUNCTION_ELEMENT_H_
#define __SUB_GRID_FUNCTION_ELEMENT_H_

#include <igatools/functions/sub_grid_function.h>
#include <igatools/functions/sub_grid_function_handler.h>
#include <igatools/functions/grid_function_element.h>

IGA_NAMESPACE_OPEN


template <int sdim,int dim,int range>
class SubGridFunctionElement
  : public GridFunctionElement<sdim,range>
{
private:
  using parent_t = GridFunctionElement<sdim,range>;
  using self_t = SubGridFunctionElement<sdim,dim,range>;

public:
  using ContainerType = const SubGridFunction<sdim,dim,range>;
  using GridElem = typename ContainerType::GridType::ElementAccessor;
  using ListIt = typename ContainerType::ListIt;

  using IndexType = typename Grid<sdim>::IndexType;

  using Value =  typename ContainerType::Value;
  template <int order>
  using Derivative = typename ContainerType::template Derivative<order>;


// using Gradient =  typename ContainerType_::Gradient;

  using Flags = grid_function_element::Flags;


  /** @name Constructors */
  ///@{
protected:
  /**
   * Default constructor. Not allowed to be used.
   */
  SubGridFunctionElement() = delete;

public:
  /**
   * Construct an accessor pointing to the element with
   * flat index @p sub_elem_index_iterator of the Function @p sub_grid_function.
   */
  SubGridFunctionElement(const std::shared_ptr<ContainerType> &sub_grid_function,
                         std::unique_ptr<GridElement<sdim>> &&sub_grid_element,
                         std::unique_ptr<GridFunctionElement<dim,range>> &&sup_grid_func_element);

  /**
   * Copy constructor. Not allowed to be used.
   */
  SubGridFunctionElement(const self_t &elem) = delete;

  /**
   * Move constructor.
   */
  SubGridFunctionElement(self_t &&elem) = default;

  /**
   * Destructor.
   */
  virtual ~SubGridFunctionElement() = default;
  ///@}


  /**
   * @name Comparison operators
   * @note In order to be meaningful, the comparison must be performed on
   * elements defined on
   * the <b>same</b> GridFunction
   * (in the sense that the pointer to the GridFunction held by the elements must
   * point to the same GridFunction object).
   */
  ///@{
  /**
   * True if the elements have the same index.
   *  @note In debug mode, it is also check they both refer to
   *  the same GridFunction. No check is done on the cache.
   */
  bool operator==(const parent_t &elem) const override;

  /**
   * True if the elements have different index.
   *  @note In debug mode, it is also check they both refer to
   *  the same GridFunction. No check is done on the cache.
   */
  bool operator!=(const parent_t &elem) const override;

  ///@}



  virtual void operator++() override;


  /**
   * @brief Move the element to the position specified by the index<tt>elem_id</tt>.
   *
   * In Debug mode an assertion will be raised
   * if the GridElement specified by <tt>elem_id</tt> has not the same property of the
   * calling GridElement.
   *
   * @warning Use this function only if you know what you are doing
   */
  void move_to(const IndexType &elem_id) override final;



  virtual void print_info(LogStream &out) const override;



  GridFunctionElement<dim,range> &
  get_sup_grid_function_element();


private:

  std::unique_ptr<GridFunctionElement<dim,range>> sup_grid_func_element_;
};


IGA_NAMESPACE_CLOSE

#endif


