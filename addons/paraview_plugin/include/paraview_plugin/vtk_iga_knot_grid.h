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
/*
class vtkPoints;
class vtkPointSet;
template<class T> class vtkSmartPointer;
//*/

IGA_NAMESPACE_OPEN

template <int dim, int codim> class Domain;

namespace paraview_plugin
{

struct VtkGridInformation;


template <class Domain>
class VtkIgaKnotGrid
{
private:

  /// Dimension.
  static const int dim = Domain::dim;

  /// Space dimension.
  static const int space_dim = Domain::space_dim;

  /**
   * Self type.
   */
  typedef VtkIgaKnotGrid Self_;

  /**
   * Self shared pointer type.
   */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /**
   * Alias for a shared pointer of a domain type.
   */
  typedef std::shared_ptr<const Domain> DomainPtr_;

  /**
   * Alias for mesh grid information shared pointer.
   */
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /**
   * Alias for vtk grid object for visualization.
   */
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /**
   * Constructor, copy and assignment operators not allowed to be used.
   */
  VtkIgaKnotGrid() = delete;
  VtkIgaKnotGrid(const VtkIgaKnotGrid &) = delete;
  VtkIgaKnotGrid(const VtkIgaKnotGrid &&) = delete;
  void operator=(const VtkIgaKnotGrid &) = delete;
  void operator=(const VtkIgaKnotGrid &&) = delete;


public:

  /**
   * Creates and returns the vtk grid for the visualization.
   */
  static VtkGridPtr_ create(const DomainPtr_ domain,
                            const GridInfoPtr_ grid_info);

private:

  /**
   * Creates and returns the vtk grid for the visualization for 1D mappings.
   */
  template<int aux_dim>
  static EnableIf<aux_dim == 1, VtkGridPtr_>
      create_grid(const DomainPtr_ domain,
                  const GridInfoPtr_ grid_info);

  /**
   * Creates and returns the vtk grid for the visualization for 2D and 3D
   *  mappings.
   */
  template<int aux_dim>
  static EnableIf<aux_dim == 2 || aux_dim == 3, VtkGridPtr_>
      create_grid(const DomainPtr_ domain,
                  const GridInfoPtr_ grid_info);

};

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_KNOT_GRID_H_
