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

#ifndef __VTK_IGA_GRID_H_
#define __VTK_IGA_GRID_H_

#include <vtkSmartPointer.h>

#include <igatools/base/config.h>

class vtkPointSet;

IGA_NAMESPACE_OPEN

class ObjectsContainer;
template <int dim, int space_dim> class IgGridFunction;

namespace paraview_plugin
{

struct VtkGridInformation;
struct VtkControlGridInformation;


template <class Domain>
class VtkIgaGrid
{
private:

  /// Dimension.
  static const int dim = Domain::dim;

  /// Space dimensions.
  static const int space_dim = Domain::space_dim;

  /// Codimension.
  static const int codim = space_dim - dim;

  /// Self type.
  typedef VtkIgaGrid Self_;

  /// Self shared pointer type.
  typedef std::shared_ptr<Self_> SelfPtr_;

  /// Alias for a shared pointer of a map function type.
  typedef std::shared_ptr<const Domain> DomainPtr_;

  /// Alias for mesh grid information shared pointer.
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /// Alias for Ig grid function of a physical grid.
  typedef IgGridFunction<Domain::dim, Domain::space_dim> IgGridFunc_;

  /**
   * Alias for vtk grid object for visualization.
   */
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /**
   * Functions container shared pointer type.
   */
  typedef std::shared_ptr<ObjectsContainer> ObjContPtr_;

  /**
   * Alias for control mesh grid information shared pointer.
   */
  using ControlGridInfoPtr_ = std::shared_ptr<VtkControlGridInformation>;

  /**
   * Constructor, copy and assignment operators not allowed to be used.
   */
  VtkIgaGrid() = delete;
  VtkIgaGrid(const Self_ &) = delete;
  VtkIgaGrid(const Self_ &&) = delete;
  Self_ &operator=(const Self_ &) = delete;
  Self_ &operator=(const Self_ &&) = delete;

  /**
   * Constructor for grids.
   */
  VtkIgaGrid(const DomainPtr_ domain,
             const Index &id,
             const GridInfoPtr_ solid_grid_info,
             const GridInfoPtr_ knot_grid_info,
             const ControlGridInfoPtr_ control_grid_info,
             const ObjContPtr_ obj_container,
             const bool is_active,
             const bool is_physical);

public:

  /**
   * TODO: to document.
   */
  static SelfPtr_ create(const DomainPtr_ domain,
                         const Index &id,
                         const GridInfoPtr_ solid_grid_info,
                         const GridInfoPtr_ knot_grid_info,
                         const ControlGridInfoPtr_ control_grid_info,
                         const ObjContPtr_ obj_container,
                         const bool is_active,
                         const bool is_physical);

  /**
   * Updates the grid information.
   */
  void update(const bool solid_updated,
              const bool knot_updated,
              const bool control_updated);

  /**
   * Returns TRUE if the object is for the physical domain.
   */
  bool is_physical() const;

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

  /**
   * Retrieves the domain name.
   */
  const std::string &get_name() const;

  /**
   * Retrieves the id.
   */
  const Index &get_id() const;

  /// Checks if the grid is active.
  bool is_active() const;

  /// Set is_active flag.
  void set_status(const bool status_flag);

  /// Checks if the grid corresponds to an ig grid funciton.
  bool is_ig_grid_func() const;


protected:

  /**
   * Mapping function for which the grids are built.
   */
  const DomainPtr_ domain_;

  /// Id of the Vtk grid.
  const Index id_;

  /**
   * Grids information for the solid grid.
   */
  GridInfoPtr_ solid_grid_info_;

  /**
   * Grids information for the knot grid.
   */
  GridInfoPtr_ knot_grid_info_;

  /**
   * Grids information for the control grid.
   */
  ControlGridInfoPtr_ control_grid_info_;

  /**
   * Container for the domain and field functions.
   */
  const ObjContPtr_ objs_container_;

  /// Vtk solid grid smart pointer.
  VtkGridPtr_ solid_grid_ = VtkGridPtr_();

  /// Vtk knot grid smart pointer.
  VtkGridPtr_ knot_grid_ = VtkGridPtr_();

  /// Vtk control grid smart pointer.
  VtkGridPtr_ control_grid_ = VtkGridPtr_();

  /// Flag for indicating if the vtk solid grid must be recomputed.
  bool recompute_solid_ = true;

  /// Flag for indicating if the vtk knot grid must be recomputed.
  bool recompute_knot_ = true;

  /// Flag for indicating if the vtk control grid must be recomputed.
  bool recompute_control_;

  /// Flag for indicating if grid is active, or not.
  bool is_active_;

  /// Flag for indicating if the generator corresponds to physical grids.
  const bool is_physical_;

  /// Flag for indicating if the generator corresponds to an ig grid func.
  const bool is_ig_grid_func_;

};

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_GRID_H_
