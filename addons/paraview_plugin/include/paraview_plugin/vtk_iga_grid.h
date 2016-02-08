//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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


/**
 * @brief This class is in charge of generating the solid, knot and
 * control polygon meshes associated to an IGA @ref Domain.
 *
 * It takes as arguments of the constructor the IGA @ref Domain,
 * the needed information for creating the VTK grids and
 * igatools objects container. Given this information,
 * it will build and return the different VTK grids (solid, knot
 * and control polygon) when required.
 *
 * @note The control polygon grids will be only created for physical domains.
 *
 * The class is aware of the possible change in the VTK grid information,
 * therefore, the different grids only will be computed (or recomputed)
 * when needed.
 *
 * For a example, if after creating a VTK grid, the type of VTK grid
 * (or the number of cells per direction) changes, the grid is recomputed.
 *
 * On the other hand, if after created, the grid is set to active,
 * when it is activated again, the grid is not recomputed.
 *
 * This considerations are done independently for the solid, not and control
 * polygon meshes.
 *
 * @todo To be defined.
 *
 * @see VtkIgaSolidGrid
 * @see VtkIgaKnotGrid
 * @see VtkIgaControlGrid
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
template <class Domain>
class VtkIgaGrid
{
private:

  /// Dimension of the @ref Domain.
  static const int dim = Domain::dim;

  /// Space dimension of the @ref Domain.
  static const int space_dim = Domain::space_dim;

  /// Codimension of the @ref Domain.
  static const int codim = space_dim - dim;

  /// Self type.
  typedef VtkIgaGrid Self_;

  /// Self shared pointer type.
  typedef std::shared_ptr<Self_> SelfPtr_;

  /// Alias for a shared pointer of the @ref Domain type.
  typedef std::shared_ptr<const Domain> DomainPtr_;

  /// Alias for VTK grid information shared pointer.
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /// Alias for VTK control grid information shared pointer.
  using ControlGridInfoPtr_ = std::shared_ptr<VtkControlGridInformation>;

  /// Alias for @ref IgGridFunction.
  typedef IgGridFunction<Domain::dim, Domain::space_dim> IgGridFunc_;

  /// Alias for VTK grid object to be created.
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /// Alias for the igatools objects container class.
  typedef std::shared_ptr<ObjectsContainer> ObjContPtr_;

  /** @name Constructors*/
  ///@{

  /**
   * @brief Constructor.
   * @param[in] domain IGA domain to be represented.
   * @param[in] id Id number associated to the VTK grid.
   * @param[in] solid_grid_info VTK grid information for the solid mesh.
   * @param[in] knot_grid_info VTK grid information for the knot mesh.
   * @param[in] control_grid_info VTK grid information for the control polygon mesh.
   * @param[in] obj_container Container storing the domains and functions
   * defined over them.
   * @param[in] Flag indicating if the VTK grid is active, or not.
   * @param[in] Flag indicating if domain is physical or parametric.
   */
  VtkIgaGrid(const DomainPtr_ domain,
             const Index &id,
             const GridInfoPtr_ solid_grid_info,
             const GridInfoPtr_ knot_grid_info,
             const ControlGridInfoPtr_ control_grid_info,
             const ObjContPtr_ obj_container,
             const bool is_active,
             const bool is_physical);

  /**
   * Default constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaGrid() = delete;

  /**
   * Copy constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaGrid(const Self_ &) = delete;

  /**
   * Move constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaGrid(const Self_ &&) = delete;
  ///@}

  /** @name Assignment operators*/
  ///@{
  /**
   * Copy assignment operator.
   * @warning Not allowed to be used.
   */
  Self_ &operator=(const Self_ &) = delete;

  /**
   * Move assignment operator.
   * @warning Not allowed to be used.
   */
  Self_ &operator=(const Self_ &&) = delete;
  ///@}

public:

  /** @name Creators*/
  ///@{
  /**
   * @brief Creates a new object of the class wrapped into a shared pointer.
   * @param[in] domain IGA domain to be represented.
   * @param[in] id Id number associated to the VTK grid.
   * @param[in] solid_grid_info VTK grid information for the solid mesh.
   * @param[in] knot_grid_info VTK grid information for the knot mesh.
   * @param[in] control_grid_info VTK grid information for the control polygon mesh.
   * @param[in] obj_container Container storing the domains and functions
   * defined over them.
   * @param[in] Flag indicating if the VTK grid is active, or not.
   * @param[in] Flag indicating if domain is physical or parametric.
   */
  static SelfPtr_ create(const DomainPtr_ domain,
                         const Index &id,
                         const GridInfoPtr_ solid_grid_info,
                         const GridInfoPtr_ knot_grid_info,
                         const ControlGridInfoPtr_ control_grid_info,
                         const ObjContPtr_ obj_container,
                         const bool is_active,
                         const bool is_physical);
  ///@}

  /**
   * @brief Updating information about which VTK grid information has changed.
   * @param[in] solid_udpated Flag indicating if the VTK grid information
   * corresponding to the solid mesh was changed.
   * @param[in] knot_udpated Flag indicating if the VTK grid information
   * corresponding to the knot mesh was changed.
   * @param[in] control_udpated Flag indicating if the VTK grid information
   * corresponding to the control polygon mesh was changed.
   */
  void update(const bool solid_updated,
              const bool knot_updated,
              const bool control_updated);

  /**
   * @brief Returns true if the domain is physical.
   * @return true if the domain is physical, false elsewhere.
   */
  bool is_physical() const;

  /**
   * @brief Computes and retrieves the VTK grid of the geometry of the domain.
   * @return Pointer to the VTK of the geometry of the domain.
   */
  VtkGridPtr_ get_solid_grid();

  /**
   * @brief Computes and retrieves the VTK grid of the knot mesh.
   * @return Pointer to the VTK of the knot mesh.
   */
  VtkGridPtr_ get_knot_grid();

  /**
   * @brief Computes and retrieves the VTK grid of the control polygon.
   * @return Pointer to the VTK of the control polygon.
   */
  VtkGridPtr_ get_control_grid();

  /**
   * @brief Retrieves the name of the grid.
   * @return Name of the grid.
   */
  const std::string &get_name() const;

  /**
   * @brief Retrieves the id number associated to the VTK grid.
   * @return Id of the grid.
   */
  const Index &get_id() const;

  /**
   * @brief Checks if the grid is active or not.
   * @return true if it is active, false elsewhere.
   */
  bool is_active() const;

  /**
   * @brief Sets and status flag.
   * @param[in] Status flag.
   */
  void set_status(const bool status_flag);

  /**
   * @brief Checks if the grid corresponds to an @ref IgGridFunction.
   * @return true, if it corresponds to an @ref IgGridFunction, false
   * elsewhere.
   */
  bool is_ig_grid_func() const;


protected:

  /// IGA domain to be represented.
  const DomainPtr_ domain_;

  /// Id associated to the VTK grid.
  const Index id_;

  /// VTK Grid information for the solid grid.
  GridInfoPtr_ solid_grid_info_;

  /// VTK Grid information for the knot grid.
  GridInfoPtr_ knot_grid_info_;

  /// VTK Grid information for the control grid.
  ControlGridInfoPtr_ control_grid_info_;

  /// Container for all the domains and function defined over them.
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

  /**
   * Flag for indicating if the generator corresponds to an
   * @ref IgGridFunction.
   */
  const bool is_ig_grid_func_;

};

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_GRID_H_
