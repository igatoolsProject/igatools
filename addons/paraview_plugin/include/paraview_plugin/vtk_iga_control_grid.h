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

#ifndef __VTK_IGA_CONTROL_GRID_H_
#define __VTK_IGA_CONTROL_GRID_H_

#include <igatools/base/config.h>

class vtkPoints;
class vtkPointSet;
template<class T> class vtkSmartPointer;

IGA_NAMESPACE_OPEN

template <int dim, int codim> class Domain;
template <int dim, int space_dim> class IgGridFunction;

namespace paraview_plugin
{

struct VtkControlGridInformation;


/**
 * @brief Helper class for creating a VTK grid representing the
 * control polygon of an IGA function.
 *
 * This class takes as argument an IGA @ref Domain, whose grid function
 * is an @ref IgGridFunction.
 *
 * This is a template class, whose argument is the @ref Domain
 * class which defines the IGA domain.
 *
 * The class extracts the control points and builds the control polygon
 * that is represented in a VTK grid, that can be structured or
 * unstructured linear.
 *
 * The control polygons can be represented for 1D, 2D and 3D IGA functions.
 *
 * The method @ref create is in charge of extracting the grid function from
 * the domain and obtaining the control points, and the methods
 * @ref create_vts_grid and @ref create_vtu_grid create the
 * structured or unstructured VTK grids.
 *
 * @note For 1D grid functions, it is not possible to create structured
 * control grids.
 *
 * @see VtkGridType
 * @see Domain
 * @see IgGridFunction
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
template <class Domain>
class VtkIgaControlGrid
{
private:
  /// Dimension of the @ref Domain.
  static const int dim = Domain::dim;

  /// Space dimension of the @ref Domain.
  static const int space_dim = Domain::space_dim;

  /// Self type of the class.
  typedef VtkIgaControlGrid Self_;

  /// Shared pointer type of the class.
  typedef std::shared_ptr<Self_> SelfPtr_;

  /// Shared pointer of the @ref Domain.
  typedef std::shared_ptr<const Domain> DomainPtr_;

  /// Shared pointer of the @ref VtkControlGridInformation.
  typedef std::shared_ptr<VtkControlGridInformation> ControlGridInfoPtr_;

  /// Shared pointer of the @ref vtkPointSet that contains the produced grid.
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /// Alias for the @ref IgGridFunction.
  typedef IgGridFunction<dim, space_dim> IgGridFun_;

  /// Shared pointer of the @ref IgGridFunction.
  typedef std::shared_ptr<const IgGridFun_> IgGridFunPtr_;

  /** @name Constructors*/
  ///@{
  /**
   * Default constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaControlGrid() = delete;

  /**
   * Copy constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaControlGrid(const VtkIgaControlGrid &) = delete;

  /**
   * Move constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaControlGrid(const VtkIgaControlGrid &&) = delete;
  ///@}

  /** @name Assignment operators*/
  ///@{
  /**
   * Copy assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkIgaControlGrid &) = delete;

  /**
   * Move assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkIgaControlGrid &&) = delete;
  ///@}


public:

  /** @name Creators*/
  ///@{
  /**
   * @brief Creates and returns the VTK grid for visualizing the control
   * grid polygon of an @ref IgGridFunction.
   *
   * This method extract the IGA function from the @ref domain and
   * extract the control points of the function.
   * It will call the method @ref create_grid_vts or @ref create_grid_vtu
   * for creating the VTK grid.
   *
   * @warning In debug mode, if the @ref domain does not corresponds to
   * a @ref IgGridFunction, an exception is thrown.
   *
   * @param[in] domain Domain defined by the IGA function whose
   * control polygon is created.
   * @param[in] grid_info Information for creating the VTK grid.
   * @return VTK grid representing the control polygon.
   */
  static VtkGridPtr_ create(const DomainPtr_ domain,
                            const ControlGridInfoPtr_ grid_info);

private:

  /**
   * @brief Creates and returns the VTK \a structured grid for visualizing the
   * control grid polygon of an @ref IgGridFunction.
   * @param[in] ig_grid_fun IGA function whose control polygon is created.
   * @param[in] points Points of the control polygon.
   * @return VTK grid representing the control polygon.
   */
  static VtkGridPtr_ create_grid_vts(const IgGridFunPtr_ ig_grid_fun,
                                     const vtkSmartPointer<vtkPoints> points);

  /**
   * @brief Creates and returns the VTK \a unstructured \a linear grid for
   * visualizing the control grid polygon of an @ref IgGridFunction.
   *
   * The unstructured grid is composed of a set of points, corresponding
   * to the control points, and lines joining them.
   *
   * @param[in] ig_grid_fun IGA function whose control polygon is created.
   * @param[in] points Points of the control polygon.
   * @return VTK grid representing the control polygon.
   */
  static VtkGridPtr_ create_grid_vtu(const IgGridFunPtr_ ig_grid_fun,
                                     const vtkSmartPointer<vtkPoints> points);
  ///@}

};

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_CONTROL_GRID_H_
