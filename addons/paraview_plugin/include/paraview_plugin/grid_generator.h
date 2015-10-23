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

#ifndef VTK_GRID_GENERATOR_H_
#define VTK_GRID_GENERATOR_H_

#include <vtkSmartPointer.h>

#include <igatools/base/config.h>


class vtkPointSet;


IGA_NAMESPACE_OPEN

template <class T, int d> class SafeSTLArray;
class FunctionsContainer;
template <int dim> class TensorSize;
template <int dim, int codim> class Domain;
template <int dim, int codim, int range, int rank> class Function;
struct VtkGridInformation;
struct VtkControlGridInformation;

template <int dim, int codim> class VtkIgaSolidGridGenerator;
template <int dim, int codim> class VtkIgaKnotGridGenerator;
template <int dim, int codim> class VtkIgaControlGridGenerator;

template <int dim, int codim>
class VtkIgaGridGenerator
{
private:
  /**
   * Space dimension.
   */
  static const int space_dim = dim + codim;

  /**
   * Self type.
   */
  typedef VtkIgaGridGenerator Self_;

  /**
   * Self shared pointer type.
   */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /**
   * Alias for a shared pointer of a map function type.
   */
  typedef std::shared_ptr<Domain<dim,codim>> DomainPtr_;

  /**
   * Alias for mesh grid information shared pointer.
   */
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /**
   * Alias for control mesh grid information shared pointer.
   */
  typedef std::shared_ptr<VtkControlGridInformation> ControlGridInfoPtr_;

  /**
   * Alias for vtk grid object for visualization.
   */
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /**
   * Functions container shared pointer type.
   */
  typedef std::shared_ptr<FunctionsContainer> FunContPtr_;

  /**
   * Vtk solid grid generator type.
   */
  typedef std::shared_ptr<VtkIgaSolidGridGenerator<dim, codim>> SolGeneratorPtr_;

  /**
   * Vtk knot grid generator type.
   */
  typedef std::shared_ptr<VtkIgaKnotGridGenerator<dim, codim>> KntGeneratorPtr_;

  /**
   * Vtk control grid generator type.
   */
  typedef std::shared_ptr<VtkIgaControlGridGenerator<dim, codim>> CtrGeneratorPtr_;

  /**
   * Constructor for grids.
   */
  VtkIgaGridGenerator(const DomainPtr_ domain,
                      const GridInfoPtr_ solid_grid_info,
                      const GridInfoPtr_ knot_grid_info,
                      const ControlGridInfoPtr_ control_grid_info,
                      const FunContPtr_ func_container,
                      const bool is_physical);

  /**
   * Constructor, copy and assignment operators not allowed to be used.
   */
  VtkIgaGridGenerator() = delete;
  VtkIgaGridGenerator(const VtkIgaGridGenerator &) = delete;
  VtkIgaGridGenerator(const VtkIgaGridGenerator &&) = delete;
  void operator=(const VtkIgaGridGenerator &) = delete;
  void operator=(const VtkIgaGridGenerator &&) = delete;

public:

  /**
   * Creates a new instance of the class and returns it wrapped into
   * a shared pointer.
   *
   * @Note: this is for physical grids.
   */
  static SelfPtr_ create_physical(const DomainPtr_ domain,
                                  const GridInfoPtr_ solid_grid_info,
                                  const GridInfoPtr_ knot_grid_info,
                                  const ControlGridInfoPtr_ control_grid_info,
                                  const FunContPtr_ func_container);

  /**
   * Creates a new instance of the class and returns it wrapped into
   * a shared pointer.
   *
   * @Note: this is for parametric grids.
   */
  static SelfPtr_ create_parametric(const DomainPtr_ domain,
                                    const GridInfoPtr_ solid_grid_info,
                                    const GridInfoPtr_ knot_grid_info,
                                    const FunContPtr_ func_container);

  /**
   * Updates the grid information for the solid, knot and control grids.
   *
   * @Note: this function is for physical grids.
   */
  void update(const bool solid_updated,
              const bool knot_updated,
              const bool control_updated);

  /**
   * Updates the grid information for the solid and knot.
   *
   * @Note: this function is for parametric grids.
   */
  void update(const bool solid_updated,
              const bool knot_updated);

  /**
   * Computes (if necessary) and returns the vtk grid the solid geometry.
   */
  VtkGridPtr_ get_solid_grid();

  /**
   * Computes (if necessary) and returns the vtk grid the knot geometry.
   */
  VtkGridPtr_ get_knot_grid();

  /**
   * Computes (if necessary) and returns the vtk grid the control geometry.
   */
  VtkGridPtr_ get_control_grid();

private:

  /**
   * Mapping function for which the grids are built.
   */
  const DomainPtr_ domain_;

  /**
   * Grids information for the solid grid.
   */
  GridInfoPtr_ solid_grid_info_;

  /**
   * Grids information for the knot grid.
   */
  GridInfoPtr_ knot_grid_info_;

  /**
   * Grids information for the knot grid.
   */
  ControlGridInfoPtr_ control_grid_info_;

  /**
   * Container for the domain and field functions.
   */
  const FunContPtr_ funcs_container_;

//    /**
//     * Number of visualization vtk cells in each direction for every Bezier
//     * element for the solid grid.
//     */
//    TensorSize<dim> n_solid_cells_;
//
//    /**
//     * Number of visualization vtk cells in each direction for every Bezier
//     * element for the knot grid.
//     */
//    TensorSize<dim> n_knot_cells_;

  /**
   * Vtk solid grid smart pointer.
   */
  VtkGridPtr_ solid_grid_;

  /**
   * Vtk knot grid smart pointer.
   */
  VtkGridPtr_ knot_grid_;

  /**
   * Vtk control grid smart pointer.
   */
  VtkGridPtr_ control_grid_;

  /**
   * Flag for indicating if the generator corresponds to physical grids,
   * or not.
   */
  bool is_physical_;

  /**
   * Flag for indicating if the vtk solid grid must be recomputed.
   */
  bool recompute_solid_;

  /**
   * Flag for indicating if the vtk knot grid must be recomputed.
   */
  bool recompute_knot_;

  /**
   * Flag for indicating if the vtk control grid must be recomputed.
   */
  bool recompute_control_;

private:
//    /**
//     * Creates a TensorSize container with the number of vtk cells in each
//     * direction.
//     */
//    static TensorSize<dim> create_n_cells (GridInfoPtr_ grid_info);

  /**
   * Computes the vtk grid the solid geometry.
   */
  void compute_solid_grid();

  /**
   * Computes the vtk grid the knot geometry.
   */
  void compute_knot_grid();

  /**
   * Computes the vtk grid the control geometry.
   */
  void compute_control_grid();

};


IGA_NAMESPACE_CLOSE

#endif // VTK_GRID_GENERATOR_H_
