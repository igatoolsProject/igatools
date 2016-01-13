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

#ifndef __VTK_IGA_SOLID_GRID_H_
#define __VTK_IGA_SOLID_GRID_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/safe_stl_vector.h>
#include <igatools/utils/tensor_size.h>
#include <paraview_plugin/vtk_iga_grid_container.h>

#include <igatools/base/objects_container.h>
#include <igatools/geometry/domain.h>
#include <igatools/functions/function.h>

class vtkPointSet;
class vtkPoints;
class vtkPointData;
class vtkCellArray;
template<class T> class vtkSmartPointer;

IGA_NAMESPACE_OPEN

class ObjectsContainer;
template <int dim> class Quadrature;
template <int dim, int codim> class Domain;

namespace paraview_plugin
{

struct VtkGridInformation;

/**
 * @brief Helper class for creating a VTK grid representing the
 * geometry of an IGA domain.
 *
 * This class takes as argument an IGA @ref Domain, that could be
 * correspond to physical or parametric domain, and represents its
 * geometry in a VTK structured or unstructured grid.
 *
 * In addition to the geometrical representation of the domain,
 * all the functions defined over the domain and contained
 * inside an igatools objects container (received as argument)
 * are visualized as point data associated to the VTK grid.
 *
 * The struct @ref PointsTopology, defined below, helps the create
 * functions in the definition the VTK grids topology.
 *
 * The grids are created by calling the public method @ref create.
 *
 * @note For 1D domains, the VTK grids are always structured.
 *
 * @author P. Antolin, 2016.
 *
 * @see Domain
 * @see Grid
 * @see GridFunction
 * @see Function
 * @see VtkGridInformation
 * @see ObjectsContainer
 *
 * @ingroup paraview_plugin
 */
template <class Domain>
class VtkIgaSolidGrid
{
private:

  /// Dimension of the @ref Domain.
  static const int dim = Domain::dim;

  /// Space dimension of the @ref Domain.
  static const int space_dim = Domain::space_dim;

  /**
   * @brief Class containing the topological information
   * for creating VTK grids representing the image of IGA domains.
   *
   * Given the grid of the domain and information of the VTK grid,
   * this class helps the class above @ref VtkIgaSolidGrid to build
   * VTK grids.
   *
   * It provides information about the number of points, cells in each
   * direction, visualization quadrature, etc.
   *
   * The main structures for defining the topology of the VTK grid are:
   *  - @ref mask_ It relate the indices of the points of the VTK
   *  cell connectivity with the indices of the quadrature points
   *  of a single IGA Bezier element.
   *  - @ref map_ For each Bezier element, it relates the index
   *  of a point of the VTK grid with a quadrature point of the Bezier
   *  element.
   *  - @ref connectivity_ Connectivity indicating the indices
   *  of the points of the VTK grid that define each cell. Only needed
   *  for VTK \a unstructured grids.
   *  - @ref quad_ Quadrature used for evaluating the points of the VTK grid.
   *  - @ref n_vis_elements_ Number of VTK point for each Bezier element
   *  in each direction.
   *
   * @see Quadrature
   * @see VtkGridInformation
   * @author P. Antolin, 2016.
   *
   * @ingroup paraview_plugin
   */
  struct PointsTopology
  {
  private:
    /// Alias for mesh grid information shared pointer.
    typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

    /// Quadrature container shared pointer type.
    typedef std::shared_ptr<Quadrature<dim>>  QuadPtr_;

    /// Connectivity container.
    typedef SafeSTLVector<SafeSTLVector<Index>>  Connectivity_;

    /// Points map container.
    typedef SafeSTLVector <SafeSTLVector <Index>> Map_;

    /// Self type;
    typedef PointsTopology Self_;

    /** @name Constructors*/
    ///@{

  public:
    /**
     * Constructor.
     * @param[in] cartesian_grid Grid for which the map and mask are computed.
     * @param[in] grid_info Information of the VTK grid.
     */
    PointsTopology(const std::shared_ptr<const Grid<dim>> cartesian_grid,
                   const GridInfoPtr_ grid_info);
  private:
    /**
     * Default constructor.
     * @warning Not allowed to be used.
     */
    PointsTopology() = delete;

    /**
     * Copy constructor.
     * @warning Not allowed to be used.
     */
    PointsTopology(const PointsTopology &) = delete;

    /**
     * Move constructor.
     * @warning Not allowed to be used.
     */
    PointsTopology(const PointsTopology &&) = delete;
    ///@}

    /** @name Assignment operators*/
    ///@{
    /**
     * Copy assignment operator.
     * @warning Not allowed to be used.
     */
    void operator=(const PointsTopology &) = delete;

    /**
     * Move assignment operator.
     * @warning Not allowed to be used.
     */
    void operator=(const PointsTopology &&) = delete;
    ///@}

  private:
    /// Connectivity for the vtk cells of a single Bezier element.
    const Connectivity_ connectivity_;

    /// Visualization quadrature.
    const QuadPtr_ quad_;

    /// Number of vtk cells per direction in each Bezier element.
    const TensorSize <dim> n_vis_elements_;

    /// Points map.
    Map_ map_;

    /// Points mask.
    SafeSTLVector <Index> mask_;

    /// Total number of points in the visualization.
    Size n_total_points_;

  public:
    /**
     * @brief Retrieves the total number of VTK visualiztion points
     * into a single IGA Bezier element.
     * @return Total number of points.
     */
    Size get_num_pts_per_bezier_elem() const;

    /**
     * @brief Retrieves the total Bezier elements into a single IGA
     * @ref Domain.
     * @return Total number of Bezier elements.
     */
    Size get_num_bezier_elems() const;

    /**
     * @brief Retrieves the total VTK cells in each Bezier element.
     * @return Total number of VTK cells into a single element.
     */
    Size get_flat_num_cells_per_bezier_elem() const;

    /**
     * @brief Retrieves the total number of points in the VTK grid.
     * @return Total number of points.
     */
    Size get_num_total_pts() const;

    /**
     * @brief Retrieves the total number of points in a single VTK cell.
     * @return Total number of points in a VTK cell.
     */
    Size get_num_pts_per_single_vtk_cell() const;

    /**
     * @brief Retrieves the mask for the points of the grid.
     * @return Mask of the grid points.
     */
    const SafeSTLVector<Index> &get_mask() const;

    /**
     * @brief Retrieves the grid points connectivity.
     * @return Grid points connectivity.
     */
    const Connectivity_ &get_connectivity() const;

    /**
     * @brief Retrieves the total number of VTK cells per Bezier element
     * along the parametric direction @p dir.
     * @param[in] dir Parametric direction.
     * @return Number of VTK cells for every Bezier element
     * along the parametric direction.
     */
    Size get_num_vtk_cells_per_bezier_elem(const Index &dir) const;

    /**
     * @brief Iterator pointing to the first position of the @ref map_.
     * @return Iterator to the beginning of the @ref map_.
     */
    Map_::const_iterator map_cbegin() const;

    /**
     * @brief Return a reference to the visualization quadrature.
     * @return Visualization quadrature.
     */
    QuadPtr_ get_quadrature() const;

  private:
    /**
     * @brief Computes and stores the points map and mask.
     * @param[in] cartesian_grid Grid for which the map and mask are computed.
     * @param[in] grid_info Information of the VTK grid.
     */
    void fill_points_map_mask(const std::shared_ptr<const Grid<dim>> cartesian_grid,
                              const GridInfoPtr_ grid_info);

    /**
     * @brief Creates the connectivity of a VTK cell for a single
     * Bezier element.
     * @param[in] grid_info Information of the VTK grid.
     * @return Connectivity of a single Bezier element.
     */
    static Connectivity_ create_element_connectivity(const GridInfoPtr_ grid_info);

    /**
     * @brief Creates the connectivity of a VTK linear cells for a single
     * Bezier element.
     * @param[in] grid_info Information of the VTK grid.
     * @return Connectivity of a single Bezier element.
     */
    static Connectivity_ create_linear_element_connectivity
    (const GridInfoPtr_ grid_info);

    /**
     * @brief Creates the connectivity of a VTK quadratic cell for a single
     * Bezier element for the 1D case.
     * @param[in] grid_info Information of the VTK grid.
     * @return Connectivity of a single Bezier element.
     */
    template <int aux_dim>
    static Connectivity_ create_quadratic_element_connectivity
    (const GridInfoPtr_ grid_info,
     typename std::enable_if_t<aux_dim == 1> * = 0);

    /**
     * @brief Creates the connectivity of a VTK quadratic cell for a single
     * Bezier element for the 2D case.
     * @param[in] grid_info Information of the VTK grid.
     * @return Connectivity of a single Bezier element.
     */
    template <int aux_dim>
    static Connectivity_ create_quadratic_element_connectivity
    (const GridInfoPtr_ grid_info,
     typename std::enable_if_t<aux_dim == 2> * = 0);

    /**
     * @brief Creates the connectivity of a VTK quadratic cell for a single
     * Bezier element for the 3D case.
     * @param[in] grid_info Information of the VTK grid.
     * @return Connectivity of a single Bezier element.
     */
    template <int aux_dim>
    static Connectivity_ create_quadratic_element_connectivity
    (const GridInfoPtr_ grid_info,
     typename std::enable_if_t<aux_dim == 3> * = 0);

    /**
     * @brief Creates and returns the visualization quadrature.
     * @param[in] grid_info Information of the VTK grid.
     * @return Visualization quadrature.
     */
    static QuadPtr_ create_visualization_quadrature(const GridInfoPtr_ grid_info);
  };

  /// Self type of the class.
  typedef VtkIgaSolidGrid Self_;

  /// Shared pointer type of the class.
  typedef std::shared_ptr<Self_> SelfPtr_;

  /// Shared pointer of the @ref Domain.
  typedef std::shared_ptr<const Domain> DomainPtr_;

  /// Shared pointer of the @ref VtkControlGridInformation.
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /// Shared pointer of the @ref vtkPointSet that contains the produced grid.
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /// Shared pointer of the quadrature type.
  typedef std::shared_ptr<Quadrature<dim>>  QuadPtr_t_;

  /// Shared pointer of the @ref ObjectsContainer.
  typedef typename std::shared_ptr<ObjectsContainer> ObjContPtr_t_;

  /**
   * Type containing all the functions with compatible dimensions for the
   * given @ref Domain.
   */
  using ValidFuncs_ = VtkIgaGridContainer::ValidFuncsForDomain<Domain>;

  /**
   * Type containing all the grid functions with compatible dimensions for
   * a the given @ref Domain assuming that corresponds to a parametric domain.
   */
  using ValidGridFuncs_ = VtkIgaGridContainer::ValidGridFuncsForDomain<Domain>;

  /** @name Constructors*/
  ///@{
  /**
   * Default constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaSolidGrid() = delete;

  /**
   * Copy constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaSolidGrid(const VtkIgaSolidGrid &) = delete;

  /**
   * Move constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaSolidGrid(const VtkIgaSolidGrid &&) = delete;
  ///@}

  /** @name Assignment operators*/
  ///@{
  /**
   * Copy assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkIgaSolidGrid &) = delete;

  /**
   * Move assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkIgaSolidGrid &&) = delete;
  ///@}


public:

  /** @name Creators*/
  ///@{
  /**
   * @brief Creates and returns the VTK grid for the visualization
   * of the domain geometry.
   *
   * In addition, the VTK grid has point data associated to it
   * that corresponds to all the function defined over the given
   * @p domain and contained in @p obj_container.
   *
   * @param[in] domain IGA domain to be represented.
   * @param[in] obj_container Container for the IGA objects.
   * @param[in] grid_info Information of the VTK grid.
   * @param[in] is_physical Boolean indicating if the domain is physical
   * (true) or parametric (false).
   */
  static VtkGridPtr_ create(const DomainPtr_ domain,
                            const ObjContPtr_t_ obj_container,
                            const GridInfoPtr_ grid_info,
                            const bool is_physical);

private:

  /**
   * @brief Creates and returns the VTK \a structured grid for the
   * visualization of the domain geometry.
   *
   * It calls the method @ref create_points for creating the points
   * of the grid, builds the structured mesh cells and creates the
   * associated point data by calling @ref create_point_data_physical or
   * @ref create_point_data_parametric.
   *
   * The point data of the VTK grid corresponds to all the functions
   * defined over the given @p domain and contained in @p obj_container.
   *
   * @param[in] domain IGA domain to be represented.
   * @param[in] obj_container Container for the IGA objects.
   * @param[in] grid_info Information of the VTK grid.
   * @param[in] is_physical Boolean indicating if the domain is physical
   * (true) or parametric (false).
   */
  static VtkGridPtr_ create_grid_vts(const DomainPtr_ domain,
                                     const ObjContPtr_t_ objs_container,
                                     const GridInfoPtr_ grid_info,
                                     const bool is_physical);

  /**
   * @brief Creates and returns the VTK \a unstructured grid for the
   * visualization of the domain geometry.
   *
   * It calls the method @ref create_points for creating the points
   * of the grid, builds the grid cells and creates the associated
   * point data by calling @ref create_point_data_physical or
   * @ref create_point_data_parametric.
   *
   * The point data of the VTK grid corresponds to all the functions
   * defined over the given @p domain and contained in @p obj_container.
   *
   * @param[in] domain IGA domain to be represented.
   * @param[in] obj_container Container for the IGA objects.
   * @param[in] grid_info Information of the VTK grid.
   * @param[in] is_physical Boolean indicating if the domain is physical
   * (true) or parametric (false).
   */
  static VtkGridPtr_ create_grid_vtu(const DomainPtr_ domain,
                                     const ObjContPtr_t_ objs_container,
                                     const GridInfoPtr_ grid_info,
                                     const bool is_physical);

  /**
   * @brief Creates and returns the points of the IGA domain needed
   * for creating the VTK grid.
   *
   * By using the information contained in @p points_top,
   * this method evaluates and stores the points of the IGA @p domain.
   *
   * @param[in] domain IGA domain to be represented.
   * @param[in] point_top Topology information for the generation of the points.
   * @return Points of the VTK mesh.
   */
  static vtkSmartPointer<vtkPoints> create_points(const DomainPtr_ domain,
                                                  const PointsTopology &points_top);

  /**
   * @brief Evaluates and associates to a given VTK grid the point
   * data corresponding to the functions defined over the given physical
   * @p domain and stored in @p objs_container.
   *
   * @param[in] domain IGA domain to be represented.
   * @param[in] obj_container Container for the IGA objects.
   * @param[in] point_top Topology information for the generation of the points.
   * @param[in,out] vtk_grid VTK grid to which the point data is added.
   */
  static void create_point_data_physical(const DomainPtr_ domain,
                                         const ObjContPtr_t_ objs_container,
                                         const PointsTopology &points_top,
                                         const VtkGridPtr_ vtk_grid);

  /**
   * @brief Evaluates and associates to a given VTK grid the point
   * data corresponding to the functions defined over the given parametric
   * @p domain and stored in @p objs_container.
   *
   * @param[in] domain IGA domain to be represented.
   * @param[in] obj_container Container for the IGA objects.
   * @param[in] point_top Topology information for the generation of the points.
   * @param[in,out] vtk_grid VTK grid to which the point data is added.
   */
  static void create_point_data_parametric(const DomainPtr_ domain,
                                           const ObjContPtr_t_ objs_container,
                                           const PointsTopology &points_top,
                                           const VtkGridPtr_ vtk_grid);
  ///@}

};

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_SOLID_GRID_H_
