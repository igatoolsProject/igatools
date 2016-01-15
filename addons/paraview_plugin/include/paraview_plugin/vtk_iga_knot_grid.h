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

#ifndef __VTK_IGA_KNOT_GRID_H_
#define __VTK_IGA_KNOT_GRID_H_

#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkIdTypeArray.h>
#include <vtkCellArray.h>

#include <igatools/base/config.h>

IGA_NAMESPACE_OPEN

template <int dim, int codim> class Domain;

namespace paraview_plugin
{

struct VtkGridInformation;

/**
 * @brief Helper class for creating a VTK grid representing the
 * knot mesh of an IGA domain.
 *
 * This class takes as argument an IGA @ref Domain, and represents
 * its knot mesh with a VTK grid. Being the knot mesh
 * the image of the knot lines of the mesh.
 *
 * This is a template class, whose argument is the @ref Domain
 * class which defines the IGA domain.
 *
 * For the case of non IGA functions, the term knot lines refers
 * to the inter-element lines.
 *
 * When we refer here to a domain, we do it in a general way: it can be
 * a parametric domain (a @ref Grid mapped with an
 * @ref grid_functions::IdentityGridFunction) or a physical domain (a
 * mapped parametric domain represented by an object of the class @ref Domain)
 * of dimensions 1, 2 and 3.
 *
 * The knot mesh is represent as a VTK unstructured grid of linear
 * or quadratic lines. It is not possible to create VTK structured grids
 * for the knot mesh.
 *
 * The grid is created by calling the public method @ref create.
 *
 * @see VtkGridType
 * @see Domain
 * @see Grid
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
template <class Domain>
class VtkIgaKnotGrid
{
private:

  /// Dimension of the @ref Domain.
  static const int dim = Domain::dim;

  /// Space dimension of the @ref Domain.
  static const int space_dim = Domain::space_dim;

  /// Self type of the class.
  typedef VtkIgaKnotGrid Self_;

  /// Shared pointer type of the class.
  typedef std::shared_ptr<Self_> SelfPtr_;

  /// Shared pointer of the @ref Domain.
  typedef std::shared_ptr<const Domain> DomainPtr_;

  /// Shared pointer of the @ref VtkGridInformation.
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /// Shared pointer of the @ref vtkPointSet that contains the produced grid.
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /** @name Constructors*/
  ///@{
  /**
   * Default constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaKnotGrid() = delete;

  /**
   * Copy constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaKnotGrid(const VtkIgaKnotGrid &) = delete;

  /**
   * Move constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaKnotGrid(const VtkIgaKnotGrid &&) = delete;
  ///@}

  /** @name Assignment operators*/
  ///@{
  /**
   * Copy assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkIgaKnotGrid &) = delete;

  /**
   * Move assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkIgaKnotGrid &&) = delete;
  ///@}

public:

  /** @name Creators*/
  ///@{
  /**
   * @brief Creates and returns the VTK grid for visualizing the knot
   * mesh of an IGA domain.
   *
   * This method calls @ref create_grid for creating the VTK grid
   * representing the knot mesh.
   *
   * @param[in] domain IGA domain to be represented.
   * @param[in] grid_info Information for creating the VTK grid.
   * @return VTK grid representing the IGA domain.
   */
  static VtkGridPtr_ create(const DomainPtr_ domain,
                            const GridInfoPtr_ grid_info);

private:

  /**
   * @brief Creates and returns the VTK grid for visualizing the knot
   * mesh of a 1D IGA domain.
   *
   * For 1D, the knot mesh is a VTK unstructured grid composed as
   * a set of points that correspond to the knot lines.
   *
   * @param[in] domain IGA domain to be represented.
   * @param[in] grid_info Information for creating the VTK grid.
   * @return VTK grid representing the 1D IGA domain.
   */
  template<int aux_dim>
  static EnableIf<aux_dim == 1, VtkGridPtr_>
  create_grid(const DomainPtr_ domain,
              const GridInfoPtr_ grid_info);

  /**
   * @brief Creates and returns the VTK grid for visualizing the knot
   * mesh of a 2D or 3D IGA domain.
   *
   * For 2D and 3D, the knot mesh is a VTK unstructured grid composed as
   * a set of edges representing the the image of the knots lines.
   *
   * @param[in] domain IGA domain to be represented.
   * @param[in] grid_info Information for creating the VTK grid.
   * @return VTK grid representing the 2D or 3D IGA domain.
   */
  template<int aux_dim>
  static EnableIf<aux_dim == 2 || aux_dim == 3, VtkGridPtr_>
  create_grid(const DomainPtr_ domain,
              const GridInfoPtr_ grid_info);
  ///@}

};

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_KNOT_GRID_H_
