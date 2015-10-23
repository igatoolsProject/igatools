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

#ifndef VTK_SOLID_GRID_GENERATOR_H_
#define VTK_SOLID_GRID_GENERATOR_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/safe_stl_vector.h>
#include <igatools/utils/tensor_size.h>

class vtkPointSet;
class vtkPoints;
class vtkPointData;
class vtkCellArray;
template<class T> class vtkSmartPointer;


IGA_NAMESPACE_OPEN

class FunctionsContainer;
template <int dim> class Quadrature;
template <int dim, int codim> class Domain;
struct VtkGridInformation;


template <int dim, int codim>
class VtkIgaSolidGridGenerator
{
private:
  /**
   * Space dimension.
   */
  static const int space_dim = dim + codim;

  /**
   * Self type.
   */
  typedef VtkIgaSolidGridGenerator Self_;

  /**
   * Self shared poitner type.
   */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /**
   * Alias for a shared pointer of a domain type.
   */
  typedef std::shared_ptr<const Domain<dim,codim>> DomainPtr_;

  /**
   * Alias for mesh grid information shared pointer.
   */
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /**
   * Alias for vtk grid object for visualization.
   */
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /**
   * Quadrature container shared pointer type.
   */
  typedef std::shared_ptr<Quadrature<dim>>  QuadPtr_t_;

  /**
   * Functions container shared pointer type.
   */
  typedef typename std::shared_ptr<FunctionsContainer> FunContPtr_t_;

  /**
   * Constructor.
   */
  VtkIgaSolidGridGenerator(const DomainPtr_ domain,
                           const GridInfoPtr_ grid_info,
                           const FunContPtr_t_ func_container);

  /**
   * Constructor, copy and assignment operators not allowed to be used.
   */
  VtkIgaSolidGridGenerator() = delete;
  VtkIgaSolidGridGenerator(const VtkIgaSolidGridGenerator &) = delete;
  VtkIgaSolidGridGenerator(const VtkIgaSolidGridGenerator &&) = delete;
  void operator=(const VtkIgaSolidGridGenerator &) = delete;
  void operator=(const VtkIgaSolidGridGenerator &&) = delete;


public:

  /**
   * Creates and returns the vtk grid for the visualization.
   */
  static VtkGridPtr_ get_grid(const DomainPtr_ domain,
                              const GridInfoPtr_ grid_info,
                              const FunContPtr_t_ func_container);

private:

  /**
   * Creates and returns the vtk grid for the visualization.
   */
  VtkGridPtr_ create_grid() const;

  /**
   * Shared pointer of the domain (i.e. the geometry).
   */
  const DomainPtr_ domain_;

  /**
   * Shared pointer of the control grid information for representing the
   * geometry.
   */
  const GridInfoPtr_ grid_info_;

  /**
   * Shared pointer of the function container.
   */
  const FunContPtr_t_ funcs_container_;

  /**
   * Number of vtk cells per direction in each Bezier element.
   */
  TensorSize <dim> n_vis_elements_;

  /**
   * Visualization quadrature.
   */
  QuadPtr_t_ quad_;

  /**
   * Points map.
   */
  SafeSTLVector <SafeSTLVector <Index>> points_map_;

  /**
   * Points mask.
   */
  SafeSTLVector <Index> points_mask_;

  /**
   * Connectivity for the vtk cells of a single Bezier element.
   */
  SafeSTLVector<SafeSTLVector<Index>> connectivity_;

  /**
   * Total number of points in the visualization.
   */
  Size n_total_points_;


  /**
   * Creates and returns the structured vtk grid for the visualization.
   */
  VtkGridPtr_ create_grid_vts() const;

  /**
   * Creates and returns the unstructured vtk grid for the visualization
   */
  VtkGridPtr_
  create_grid_vtu() const;

  /**
   * Creates the points for the visualization grids.
   * (Both, structured and unstructured).
   */
  vtkSmartPointer<vtkPoints> create_points() const;

  /**
   * Initializes the points map and mask, and the total number of points in
   * the visualization.
   */
  void init_points_info();

  /**
   * Creates the connectivity of the VTK linear cells for a single
   * Bezier element.
   */
  void create_linear_element_connectivity();

  /**
   * Creates the connectivity of the VTK quadratic cells for a single
   * Bezier element. For the 1D case.
   */
  template <int aux_dim>
  void create_quadratic_element_connectivity
  (typename std::enable_if_t<aux_dim == 1> * = 0);

  /**
   * Creates the connectivity of the VTK quadratic cells for a single
   * Bezier element. For the 2D case.
   */
  template <int aux_dim>
  void create_quadratic_element_connectivity
  (typename std::enable_if_t<aux_dim == 2> * = 0);

  /**
   * Creates the connectivity of the VTK quadratic cells for a single
   * Bezier element. For the 3D case.
   */
  template <int aux_dim>
  void create_quadratic_element_connectivity
  (typename std::enable_if_t<aux_dim == 3> * = 0);

  /**
   * Creates the quadrature needed for the visualization.
   */
  void create_visualization_quadrature();

  /**
   * Creates the point data associated to the mapping.
   */
  template <int range, int rank>
  void create_point_data(vtkPointData *const point_data) const;

  /**
   * Auxiliar method for creating the point data associated to the mapping
   * for dim = 1 and codim = 0.
   */
  template <int aux_dim, int aux_codim>
  void
  create_point_data_dim_codim(vtkPointData *const data,
                              typename std::enable_if_t<(aux_dim == 1 && aux_codim == 0)>* = 0) const;

  /**
   * Auxiliar method for creating the point data associated to the mapping
   * for the cases that are not dim = 1 and codim = 0.
   */
  template <int aux_dim, int aux_codim>
  void
  create_point_data_dim_codim(vtkPointData *const data,
                              typename std::enable_if_t<!(aux_dim == 1 && aux_codim == 0)>* = 0) const;


  /**
   * Copies data from a igatools tensor to a vtk tuple like structure.
   */
  void tensor_to_tuple(const Tdouble t, Real *const tuple, int &pos) const;

  /**
   * Copies data from a igatools tensor to a vtk tuple like structure.
   */
  template <class Tensor>
  void tensor_to_tuple(const Tensor &t, Real *const tuple, int &pos) const;

  /**
   * Copies data from a igatools tensor to a vtk tuple like structure.
   */
  template <class Tensor>
  void tensor_to_tuple(const Tensor &t, Real *const tuple) const;

};

IGA_NAMESPACE_CLOSE

#endif // VTK_SOLID_GRID_GENERATOR_H_
