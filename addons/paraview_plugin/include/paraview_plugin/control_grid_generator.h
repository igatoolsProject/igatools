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

#ifndef __VTK_CONTROL_GRID_GENERATOR_H_
#define __VTK_CONTROL_GRID_GENERATOR_H_

#include <igatools/base/config.h>

class vtkPoints;
class vtkPointSet;
template<class T> class vtkSmartPointer;

IGA_NAMESPACE_OPEN

template <int dim, int space_dim> class IgGridFunction;
template <int dim, int codim> class Domain;
struct VtkControlGridInformation;


template <class Domain>
class VtkIgaControlGridGenerator
{
private:
  /// Dimension.
  static const int dim = Domain::dim;

  /// Space dimension.
  static const int space_dim = Domain::space_dim;

  /**
   * Self type.
   */
  typedef VtkIgaControlGridGenerator Self_;

  /**
   * Self shared poitner type.
   */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /**
   * Alias for a shared pointer of a map function type.
   */
  typedef std::shared_ptr<const Domain> DomainPtr_;

  /**
   * Alias for mesh grid information shared pointer.
   */
  typedef std::shared_ptr<VtkControlGridInformation> ControlGridInfoPtr_;

  /**
   * Alias for vtk grid object for visualization.
   */
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /**
   * Alias for a ig function type.
   */
  typedef IgGridFunction<dim, space_dim> IgGridFun_;

  /**
   * Alias for a shared pointer of a ig function type.
   */
  typedef std::shared_ptr<const IgGridFun_> IgGridFunPtr_;

  /**
   * Constructor, copy and assignment operators not allowed to be used.
   */
  VtkIgaControlGridGenerator() = delete;
  VtkIgaControlGridGenerator(const VtkIgaControlGridGenerator &) = delete;
  VtkIgaControlGridGenerator(const VtkIgaControlGridGenerator &&) = delete;
  void operator=(const VtkIgaControlGridGenerator &) = delete;
  void operator=(const VtkIgaControlGridGenerator &&) = delete;


public:

  /**
   * Creates and returns the vtk grid for the visualization.
   */
  static VtkGridPtr_ create(const DomainPtr_ domain,
                            const ControlGridInfoPtr_ grid_info);

private:

  /**
   * Creates and returns the vtk structured grid for the visualization.
   */
  static VtkGridPtr_ create_grid_vts(const IgGridFunPtr_ ig_grid_fun,
                                     const vtkSmartPointer<vtkPoints> points);

  /**
   * Creates and returns the vtk unstructured grid for the visualization.
   */
  static VtkGridPtr_ create_grid_vtu(const IgGridFunPtr_ ig_grid_fun,
                                     const vtkSmartPointer<vtkPoints> points);

};

IGA_NAMESPACE_CLOSE

#endif // __VTK_CONTROL_GRID_GENERATOR_H_
