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

#ifndef __VTK_IGA_GRID_INFORMATION_H_
#define __VTK_IGA_GRID_INFORMATION_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_size.h>

#include "vtk_iga_types.h"

IGA_NAMESPACE_OPEN

namespace paraview_plugin
{

/**
 * @brief @p Struct containing information for building a VTK grid for
 * IGA solid and knot grids.
 *
 * This struct contains the basic options related to the VTK grid
 * for building @ref VtkIgaSolidGrid and @ref VtkIgaKnotGrid.
 *
 * The information is:
 *  - @ref grid_type_ that is a @ref VtkGridType of the grid.
 *  - @ref cells_per_element_ the number of VTK cells in each
 *    parametric direction for each Bezier element of an IGA domain.
 *    The number of cells in each direction must be greather than 0 (on
 *    the contrary, an exception will be thrown).
 *
 * The number of VTK cells in each direction is defined for the 3 physical
 * directions. For 1D and 2D domains, only the first or the two first
 * values will be used.
 *
 * The struct provides methods for retrieving the information.
 *
 * @see VtkGridType
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
struct VtkGridInformation
{
private:
  // Container type for the number of VTK cells per Bezier element.
  typedef SafeSTLArray<int, 3> NumCellsContainer_;

  /// Self type of the class.
  typedef VtkGridInformation Self_;

  /// Shared pointer type of the class.
  typedef std::shared_ptr<Self_> SelfPtr_;

  /** @name Constructors*/
  ///@{
  /**
   * Constructor.
   * @param[in] num_cells Number of VTK cells in each direction
   * for each Bezier element.
   * @param[in] grid_type Type of the VTK grid.
   */
  VtkGridInformation(const NumCellsContainer_ &num_cells,
                     const VtkGridType &grid_type);

  /**
   * Default constructor.
   * @warning Not allowed to be used.
   */
  VtkGridInformation() = delete;

  /**
   * Copy constructor.
   * @warning Not allowed to be used.
   */
  VtkGridInformation(const VtkGridInformation &grid_info) = delete;

  /**
   * Move constructor.
   * @warning Not allowed to be used.
   */
  VtkGridInformation(const VtkGridInformation &&) = delete;
  ///@}

  /** @name Assignment operators*/
  ///@{

  /**
   * Copy assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkGridInformation &grid_info) = delete;

  /**
   * Move assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkGridInformation &&) = delete;
  ///@}

public:

  /** @name Creators*/
  ///@{

  /**
   * Creates and returns a shared pointer with a new instance of the class.
   *
   * @param[in] num_cells Number of VTK cells in each direction
   * for each Bezier element.
   * @param[in] grid_type Type of the VTK grid.
   * @return New instance of the class wrapped inside a shared pointer.
   */
  static SelfPtr_ create(const NumCellsContainer_ &num_cells,
                         const VtkGridType &grid_type);

  /**
   * Creates and returns a shared pointer with a default new instance of
   * the class.
   *
   * The number of VTK cells in each direction will be set to 1,
   * and the @ref VtkGridType will be set to @p none.
   *
   * @return New instance of the class wrapped inside a shared pointer.
   */
  static SelfPtr_ create_void();
  ///@}

  /**
   * @brief Updates the information.
   *
   * It receives a @VtkGridInformation object for updating the
   * information of the class.
   * If the information of the argument is equal to the information
   * of the class, it is not updated.
   *
   * @param[in] @VtkGridInformation object for updating the
   * information of the class.
   * @return True if the information was updated, false elsewhere.
   */
  bool update(SelfPtr_ grid_info);

  /**
   * @brief Retrieves the grid type.
   * @return VTK grid type.
   */
  const VtkGridType &get_grid_type() const;

  /**
   * @brief Returns true or false if the VTK grid is structured, or not.
   * @return True if the VTK grid type is structured, false elsewhere.
   */
  bool is_structured() const;

  /**
   * @brief Returns true or false if the VTK grid is made of quadratic
   * cells, or not.
   * @return True if the VTK grid type is quadratic, false elsewhere.
   */
  bool is_quadratic() const;

  /**
   * TODO: to document.
   * Returns the number of cells per element.
   */
  /**
   * @brief Returns the number of VTK cells per Bezier element
   * for the given dimension @p dim.
   * @tparam dim Dimension of the returning array.
   * @return Number of cells per element.
   */
  template <int dim>
  SafeSTLArray <Size, dim> get_num_cells_per_element() const;

  /**
   * @brief Prints the information of the class for debugging purposes.
   * @param[in] out @ref LogStream for printing the information.
   */
  void print_info(LogStream &out) const;

private:

  /// VTK grid type.
  VtkGridType grid_type_;

  /// Number of cells in each parametric direction in each Bezier element.
  NumCellsContainer_ cells_per_element_;
};



/**
 * @brief @p Struct containing information for building a VTK grid
 * for IGA control point grids.
 *
 * The information is stored in @ref grid_type_ that is an object.
 * of the type @ref VtkGridType.
 *
 * Control point grids can only be represented by means of VTK
 * structured grids, or VTK unstructured linear grids.
 *
 * The struct provides methods for retrieving the information.
 *
 * @see VtkGridType
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
struct VtkControlGridInformation
{
private:

  /// Self type of the class.
  typedef VtkControlGridInformation Self_;

  /// Shared pointer type of the class.
  typedef std::shared_ptr<Self_> SelfPtr_;

  /** @name Constructors*/
  ///@{
  /**
   * Constructor.
   * @param[in] structured For indicating if the VTK grid is structure,
   * of unstructured made of linear cells.
   */
  VtkControlGridInformation(const bool structured);

public:
  /**
   * Default constructor.
   * It considers the @ref VtkGridType as @p structured.
   */
  VtkControlGridInformation() : VtkControlGridInformation(true) {};

private:
  /**
   * Copy constructor.
   * @warning Not allowed to be used.
   */
  VtkControlGridInformation(const VtkControlGridInformation &) = delete;

  /**
   * Move constructor.
   * @warning Not allowed to be used.
   */
  VtkControlGridInformation(const VtkControlGridInformation &&) = delete;
  ///@}

  /** @name Assignment operators*/
  ///@{
  /**
   * Copy assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkControlGridInformation &) = delete;

  /**
   * Move assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkControlGridInformation &&) = delete;
  ///@}


public:
  /** @name Creators*/
  ///@{

  /**
   * Creates and returns a shared pointer with a new instance of the class.
   *
   * @param[in] structured For indicating if the VTK grid is structure,
   * of unstructured made of linear cells.
   * @return New instance of the class wrapped inside a shared pointer.
   */
  static SelfPtr_ create(const bool structured);

  /**
   * Creates and returns a shared pointer with a default new instance of
   * the class.
   *
   * The @ref VtkGridType is set to @p structured.
   *
   * @return New instance of the class wrapped inside a shared pointer.
   */
  static SelfPtr_ create_void();
  ///@}

  /**
   * @brief Updates the information.
   *
   * It receives a @VtkGridInformation object for updating the
   * information of the class.
   * If the information of the argument is equal to the information
   * of the class, it is not updated.
   *
   * @param[in] @VtkGridInformation object for updating the
   * information of the class.
   * @return True if the information was updated, false elsewhere.
   */
  bool update(SelfPtr_ grid_info);

  /**
   * @brief Returns true or false if the VTK grid is structured, or not.
   * @return True if the VTK grid type is structured, false elsewhere.
   */
  bool is_structured() const;

  /**
   * @brief Retrieves the grid type.
   * @return VTK grid type.
   */
  const VtkGridType &get_grid_type() const;

  /**
   * @brief Prints the information of the class for debugging purposes.
   * @param[in] out @ref LogStream for printing the information.
   */
  void print_info(LogStream &out) const;

private:

  /// Vtk grid type.
  VtkGridType grid_type_;

};

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_GRID_INFORMATION_H_
