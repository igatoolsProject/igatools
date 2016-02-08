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

#ifndef __VTK_IGA_GRID_CONTAINER_H_
#define __VTK_IGA_GRID_CONTAINER_H_

#include <igatools/base/config.h>

#include <paraview_plugin/vtk_iga_types.h>
#include <igatools/base/objects_container.h>
#include <igatools/geometry/grid.h>
#include <igatools/functions/grid_function.h>
#include <igatools/geometry/domain.h>
#include <igatools/functions/function.h>

class vtkMultiBlockDataSet;

IGA_NAMESPACE_OPEN

class ObjectsContainer;

namespace paraview_plugin
{

template <class Domain> class VtkIgaGrid;
struct VtkGridInformation;
struct VtkControlGridInformation;


/**
 * @brief Creates the VTK grids for all the IGA domains contained
 * in the input file.
 *
 * The main purpose of this class is to parse an objects container defined
 * in an input file and create, when required, the VTK grids of the
 * solid, knot and control polygon meshes for all the physical and
 * parametric domains in the container.
 *
 * All the created will be inserted into a single VTK multiblock data set.
 * The structure of the multiblock is:
 *  - Multiblock
 *      - Physical domains
 *          - Solid meshes
 *          - Knot meshes
 *          - Control polygon meshes
 *      - Parametric domains
 *          - Solid meshes
 *          - Knot meshes
 *
 * All the physical domains will be created by means of a group
 * of @ref VtkGridType stored in @ref phys_grids_. For the parametric
 * domains, @ref parm_grids_ will create the grids.
 *
 * It takes as input the input file, the VTK information about
 * the number of cells per direction and the @ref VtkGridType
 * for each type mesh.
 *
 * This class also provides methods for communicating with ParaView.
 * Thus, the different domains (physical or parametric) can be activated
 * or de-activated (they will be visualized or not).
 *
 * It also provide methods for querying how many domains exists and
 * querying their names.
 *
 * Currently, the plugin is able to visualize domains whose parametric
 * and space dimensions are 1, 2 or 3 (that hereinafter we call them
 * \a valid \a domains).
 * 0D entities are not supported yet.
 * ParaView is not able to visualize 4D (or higher) geometrical objects.
 *
 * @note The \a valid domains and functions that can be visualized are
 * a subset of the instantiated dimensions in igatools.
 *
 * Therefore, warnings will be thrown if \a invalid \a domains are found.
 *
 * In addition to the representation of the domains, this class also takes
 * care of the functions defined over those domains to be visualized also
 * in ParaView as data associated to the grids.
 *
 * Once the class has been constructed, the VTK grid information can
 * be modified. In such case, only those grids that need to be recomputed
 * will be created again.
 *
 * On the other hand, it is not possible to update the objects container
 * (by means of a new input file). In that case, a new class must be
 * created.
 *
 * The input file can be written in ascii format as a XML human readable
 * objects container file (see @subpage in_out for further details);
 * or in binary format as the result of the serialization of an objects
 * container.
 *
 * @warning Cases have been found in which trying to de-serialize
 * an objects container from a binary file (that has not been written
 * properly) the de-serialize operation stacks, instead of throwing
 * an error that could be managed.
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
class VtkIgaGridContainer
{
private:

  /// Self type of the class.
  typedef VtkIgaGridContainer Self_;

  /// Shared pointer type of the class.
  typedef std::shared_ptr<Self_> SelfPtr_;

  /// Container for the number of VTK cells by Bezier element in each direction.
  typedef SafeSTLArray<int, 3> NumCells_;

  /**
   * @brief Struct for determining if a type @p T has a wrong parametric
   * dimension.
   *
   * Valid dimensions are 1, 2 or 3.
   */
  template <class T>
  struct IsInValidDim :
    boost::mpl::or_<
  boost::mpl::bool_<(T::dim < 1)>,
  boost::mpl::bool_<(T::dim > 3)>> {};

  /**
   * @brief Struct for determining if a type @p T has a wrong space
   * dimension.
   *
   * Valid space dimensions are 1, 2 or 3.
   */
  template <class T>
  struct IsInValidSpaceDim :
    boost::mpl::or_<
  boost::mpl::bool_<(T::space_dim < 1)>,
  boost::mpl::bool_<(T::space_dim > 3)>> {};

  /**
   * @brief Struct for determining if a type @ref Domain has invalid
   * dimensions.
   *
   * Valid domains are those which have their parametric and space
   * dimensions valid according to the definition of @ref IsInValidDim and
   * @ref IsInValidSpaceDim.
   */
  template <class Domain>
  struct IsInValidDomain :
    boost::mpl::or_<IsInValidDim<Domain>, IsInValidSpaceDim<Domain>> {};

  /**
   * @brief Struct for determining if a type @ref Function has invalid
   * dimensions.
   *
   * Valid functions are those which have their parametric and space
   * dimensions valid according to the definition of @ref IsInValidDim and
   * @ref IsInValidSpaceDim.
   */
  template <class Function>
  struct IsInValidFunction :
    boost::mpl::or_<IsInValidDim<Function>, IsInValidSpaceDim<Function>> {};

  /**
   * @brief Struct for determining if the dimensions of a @ref Function
   * match with the @ref Domain dimensions.
   *
   * The parametric and space dimensions must be equal.
   */
  template <class Function, class Domain>
  struct IsInValidFunctionForDomain_ :
    boost::mpl::or_<
  boost::mpl::bool_<(Function::dim != Domain::dim)>,
  boost::mpl::bool_<(Function::space_dim != Domain::space_dim)>> {};

  /**
   * @brief Struct for determining if the dimensions of a @ref GridFunction
   * match with the @ref Domain (a parametric domain) dimensions
   *
   * The parametric dimension must be equal for both, and, in addition,
   * the codimention of the @ref Domain must be 0 (in order to correspond
   * to a parametric domain).
   */
  template <class GridFunction, class Domain>
  struct IsInValidGridFunctionForDomain_ :
    boost::mpl::or_<
  boost::mpl::bool_<(GridFunction::dim != Domain::dim)>,
  boost::mpl::bool_<(Domain::space_dim != Domain::dim)>> {};

  /**
   * Alias for a <tt>boost::mpl vector</tt> containing all the valid
   * @ref Domain classes for the plugin.
   */
  using ValidDomains_ = boost::mpl::remove_if<
                        InstantiatedTypes::Domains,
                        boost::mpl::lambda< IsInValidDomain< boost::mpl::_1 > >::type
                        >::type;

  /**
   * Alias for a <tt>boost::mpl vector</tt> containing all the valid
   * @ref Grid classes for the plugin.
   */
  using ValidGrids_ = boost::mpl::remove_if<
                      InstantiatedTypes::Grids,
                      boost::mpl::lambda< IsInValidDim< boost::mpl::_1 > >::type
                      >::type;

  /**
   * Alias for a <tt>boost::mpl vector</tt> containing all the valid
   * @ref GridFunction classes for the plugin.
   */
  using ValidGridFuncs_ = typename boost::mpl::remove_if<
                          InstantiatedTypes::GridFunctions,
                          typename boost::mpl::lambda< IsInValidDim< boost::mpl::_1 > >::type
                          >::type;
  /**
   * Alias for a <tt>boost::mpl vector</tt> containing all the valid
   * @ref Function classes for the plugin.
   */
  using ValidFunctions_ = typename boost::mpl::remove_if<
                          InstantiatedTypes::Functions,
                          typename boost::mpl::lambda< IsInValidDomain< boost::mpl::_1 > >::type
                          >::type;

  /**
   * Alias for a <tt>boost::mpl vector</tt> containing all the \a invalid
   * @ref Domain classes for the plugin.
   */
  using InvalidDomains_ = boost::mpl::remove_if<
                          InstantiatedTypes::Domains,
                          boost::mpl::lambda< boost::mpl::not_<IsInValidDomain< boost::mpl::_1 >> >::type
                          >::type;

  /**
   * Alias for a <tt>boost::mpl vector</tt> containing all the \a invalid
   * @ref Grid classes for the plugin.
   */
  using InvalidGrids_ = boost::mpl::remove_if<
                        InstantiatedTypes::Grids,
                        boost::mpl::lambda< boost::mpl::not_<IsInValidDim< boost::mpl::_1 >> >::type
                        >::type;

  /**
   * Alias for a <tt>boost::mpl vector</tt> containing all the \a invalid
   * @ref GridFunction classes for the plugin.
   */
  using InvalidGridFuncs_ = typename boost::mpl::remove_if<
                            InstantiatedTypes::GridFunctions,
                            boost::mpl::lambda< boost::mpl::not_<IsInValidDim< boost::mpl::_1 >> >::type
                            >::type;

  /**
   * Alias for a <tt>boost::mpl vector</tt> containing all the \a invalid
   * @ref Function classes for the plugin.
   */
  using InvalidFunctions_ = typename boost::mpl::remove_if<
                            InstantiatedTypes::Functions,
                            boost::mpl::lambda< boost::mpl::not_<IsInValidDomain< boost::mpl::_1 >> >::type
                            >::type;


  /**
   * @brief This functor transform a given <tt>boost::mpl vector T</tt>
   * into a <tt>boost::fusion vector</tt>.
   *
   * The new <tt>boost::fusion vector</tt> created contains a shared pointer
   * of each type contained in the vector @p T.
   */
  template< class T >
  struct AsFusionVectorSharedPtr
  {
    typedef typename boost::fusion::result_of::as_vector<
    typename boost::mpl::transform<T, std::shared_ptr<boost::mpl::_1>>::type>::type type;
  };

  /// <tt>boost::fusion vector</tt> containing the types of all the valid grids.
  using GridPtrs_ = AsFusionVectorSharedPtr<ValidGrids_>::type;

  /// <tt>boost::fusion vector</tt> containing the types of all the valid grids functions.
  using GridFuncPtrs_ = AsFusionVectorSharedPtr<ValidGridFuncs_>::type;

  /// <tt>boost::fusion vector</tt> containing the types of all the valid domains.
  using DomainPtrs_ = AsFusionVectorSharedPtr<ValidDomains_>::type;

  /// <tt>boost::fusion vector</tt> containing the types of all the valid functions.
  using FunctionPtrs_ = AsFusionVectorSharedPtr<ValidFunctions_>::type;

  /// <tt>boost::fusion vector</tt> containing the types of all the \a invalid grids.
  using InvalidGridPtrs_ = AsFusionVectorSharedPtr<InvalidGrids_>::type;

  /// <tt>boost::fusion vector</tt> containing the types of all the \a invalid grids functions.
  using InvalidGridFuncPtrs_ = AsFusionVectorSharedPtr<InvalidGridFuncs_>::type;

  /// <tt>boost::fusion vector</tt> containing the types of all the \a invalid domains.
  using InvalidDomainPtrs_ = AsFusionVectorSharedPtr<InvalidDomains_>::type;

  /// <tt>boost::fusion vector</tt> containing the types of all the \a invalid functions.
  using InvalidFunctionPtrs_ = AsFusionVectorSharedPtr<InvalidFunctions_>::type;

public:
  /**
   * Valid functions for a given domain.
   */
  template <class Domain>
  using ValidFuncsForDomain = typename boost::fusion::result_of::as_vector<
                              typename boost::mpl::transform<
                              typename boost::mpl::remove_if<
                              InstantiatedTypes::Functions,
                              typename boost::mpl::lambda< IsInValidFunctionForDomain_< boost::mpl::_1, Domain > >::type>::type,
                              std::shared_ptr<boost::mpl::_1>>::type>::type;

  /**
   * Valid grid functions for a given domain.
   */
  template <class Domain>
  using ValidGridFuncsForDomain = typename boost::fusion::result_of::as_vector<
                                  typename boost::mpl::transform<
                                  typename boost::mpl::remove_if<
                                  InstantiatedTypes::GridFunctions,
                                  typename boost::mpl::lambda< IsInValidGridFunctionForDomain_< boost::mpl::_1, Domain > >::type>::type,
                                  std::shared_ptr<boost::mpl::_1>>::type>::type;

private:

  /**
   * @brief This functor transform a given <tt>boost::mpl vector T</tt>
   * into a <tt>boost::fusion map</tt>.
   *
   * The new <tt>boost::fusion map</tt> created contains, for each type
   * of @p T, a pair composed by the type @p T itself
   * and a vector of shared pointers @ref VtkIgaGrid for the type @p T.
   */
  template< class T >
  struct AsFusionMapTypeVtkIgaGridVector
  {
  private:
    /// Alias of a shared pointer of the VTK IGA grid type for a given @ref Domain type.
    template <class Domain>
    using Pair_ = boost::fusion::pair<Domain,
          SafeSTLVector<std::shared_ptr<VtkIgaGrid<Domain>>>>;

  public:
    typedef typename boost::fusion::result_of::as_map<
    typename boost::mpl::transform<T, Pair_<boost::mpl::_1>>::type>::type type;
  };

  /**
   * Basic container for VTK IGA grids.
   */
  using VtkIgaGridsContainer_ = AsFusionMapTypeVtkIgaGridVector<ValidDomains_>::type;

  /// Grid information shared pointer type.
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /// Objects container shared pointer type.
  typedef std::shared_ptr<ObjectsContainer> ObjContPtr_;

  /// Control grid information shared pointer type.
  typedef std::shared_ptr<VtkControlGridInformation> ControlGridInfoPtr_;

  /** @name Constructors*/
  ///@{
  /**
   * Constructor constructor.
   * @param[in] objs_container Container of all the domains to be
   * represented and it associated functions.
   * @param[in] phys_solid_info VTK grid information for the solid meshes
   * of physical domains.
   * @param[in] phys_knot_info VTK grid information for the knot meshes
   * of physical domains.
   * @param[in] phys_control_info VTK grid information for the control
   * polygon meshes of physical domains.
   * @param[in] parm_solid_info VTK grid information for the solid meshes
   * of parametric domains.
   * @param[in] parm_knot_info VTK grid information for the knot meshes
   * of parametric domains.
   */
  VtkIgaGridContainer(const ObjContPtr_ objs_container,
                      const GridInfoPtr_ phys_solid_info,
                      const GridInfoPtr_ phys_knot_info,
                      const ControlGridInfoPtr_ phys_control_info,
                      const GridInfoPtr_ parm_solid_info,
                      const GridInfoPtr_ parm_knot_info);

  /**
   * Void constructor.
   */
  VtkIgaGridContainer();

  /**
   * Copy constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaGridContainer(const VtkIgaGridContainer &) = delete;

  /**
   * Move constructor.
   * @warning Not allowed to be used.
   */
  VtkIgaGridContainer(const VtkIgaGridContainer &&) = delete;
  ///@}

  /** @name Assignment operators*/
  ///@{

  /**
   * Copy assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkIgaGridContainer &) = delete;

  /**
   * Move assignment operator.
   * @warning Not allowed to be used.
   */
  void operator=(const VtkIgaGridContainer &&) = delete;
  ///@}

public:

  /** @name Creators*/
  ///@{
  /**
   * @brief This methods creates a new instance of the class and returns
   * it wrapped into a shared pointer.
   *
   * It parses an objects container from the given file (by calling
   * @ref parse_objects_container) and takes the VTK grid options passed
   * from @ref IgatoolsParaViewReader for creating the
   * @ref VtkIgaGridInformation and @ref VtkIgaControlGridInformation for
   * all the grid types.
   *
   * @param[in] file_name Input file containing the objects container.
   * @param[in] n_cells_phs_sol Number of VTK cells in each parametric
   * direction for each Bezier element of the physical domain solid meshes.
   * @param[in] n_cells_phs_knt Number of VTK cells in each parametric
   * direction for each Bezier element of the physical domain knot meshes.
   * @param[in] n_cells_phs_ctr Number of VTK cells in each parametric
   * direction for each Bezier element of the physical domain control
   * polygon meshes.
   * @param[in] n_cells_prm_sol Number of VTK cells in each parametric
   * direction for each Bezier element of the parametric domain solid meshes.
   * @param[in] n_cells_phs_knt Number of VTK cells in each parametric
   * direction for each Bezier element of the parametric domain knot meshes.
   *
   * @return New instance of the class wrapped into a shared pointer.
   */
  static SelfPtr_ create(const std::string &file_name,
                         const NumCells_   &n_cells_phs_sol,
                         const VtkGridType &grid_type_phs_sol,
                         const NumCells_   &n_cells_phs_knt,
                         const VtkGridType &grid_type_phs_knt,
                         const VtkGridType &grid_type_phs_ctr,
                         const NumCells_   &n_cells_prm_sol,
                         const VtkGridType &grid_type_prm_sol,
                         const NumCells_   &n_cells_prm_knt,
                         const VtkGridType &grid_type_prm_knt);

  /**
   * @brief This methods creates a new empty instance of the class and
   * returns it wrapped into a shared pointer.
   *
   * All the options are set to default, and the objects container is
   * empty, so, it will not be able to create any geometry.
   * @return New instance of the class wrapped into a shared pointer.
   */
  static SelfPtr_ create_void();

  ///@}

private:

  /// Objects container.
  ObjContPtr_ objs_container_;

  /// Grids information for the physical solid mesh.
  const GridInfoPtr_ phys_solid_info_;

  /// Grids information for the physical knot mesh.
  const GridInfoPtr_ phys_knot_info_;

  /// Grids information for the physical control mesh.
  const ControlGridInfoPtr_ phys_control_info_;

  /// Grids information for the parametric solid mesh.
  const GridInfoPtr_ parm_solid_info_;

  /// Grids information for the parametric knot mesh.
  const GridInfoPtr_ parm_knot_info_;

  /// Collection of physical VTK IGA grids.
  VtkIgaGridsContainer_ phys_grids_;

  /// Collection of parametric VTK IGA grid.
  VtkIgaGridsContainer_ parm_grids_;


public:
  /**
   * @brief Updates the VTK grid information for all grid types.
   *
   * Once all the @ref VtkGridInformation and @ref VtkControlGridInformation
   * have been updated, all the @ref VtkIgaGrid (all the ones
   * contained in @ref phys_grids_ and @ref parm_grids_) are
   * also updated.
   *
   * @param[in] n_cells_phs_sol Number of VTK cells in each parametric
   * direction for each Bezier element of the physical domain solid meshes.
   * @param[in] n_cells_phs_knt Number of VTK cells in each parametric
   * direction for each Bezier element of the physical domain knot meshes.
   * @param[in] n_cells_phs_ctr Number of VTK cells in each parametric
   * direction for each Bezier element of the physical domain control
   * polygon meshes.
   * @param[in] n_cells_prm_sol Number of VTK cells in each parametric
   * direction for each Bezier element of the parametric domain solid meshes.
   * @param[in] n_cells_phs_knt Number of VTK cells in each parametric
   * direction for each Bezier element of the parametric domain knot meshes.
   */
  void update(const NumCells_   &n_cells_phs_sol,
              const VtkGridType &grid_type_phs_sol,
              const NumCells_   &n_cells_phs_knt,
              const VtkGridType &grid_type_phs_knt,
              const VtkGridType &grid_type_phs_ctr,
              const NumCells_   &n_cells_prm_sol,
              const VtkGridType &grid_type_prm_sol,
              const NumCells_   &n_cells_prm_knt,
              const VtkGridType &grid_type_prm_knt);

  /**
   * @brief Computes the required VTK grid and fills them into the
   * multiblock data set.
   *
   * Given indicating which grid types are active, this method creates
   * the required VTK grids and inserts them into the multiblock structure.
   *
   * @param[in] phys_mesh Flag indicating if the physical mesh is active.
   * @param[in] sol_phys_mesh Flag indicating if the solid physical mesh is active.
   * @param[in] knot_phys_mesh Flag indicating if the knot physical mesh is active.
   * @param[in] ctr_phys_mesh Flag indicating if the control polygon physical mesh is active.
   * @param[in] prm_mesh Flag indicating if the parametric mesh is active.
   * @param[in] sol_prm_mesh Flag indicating if the solid parametric mesh is active.
   * @param[in] knot_prm_mesh Flag indicating if the knot parametric mesh is active.
   * @param[out] mb VTK multiblock data set to be filled with the corresponding grids.
   */
  void create_multiblock_grid(const bool phys_mesh,
                              const bool sol_phys_mesh,
                              const bool knt_phys_mesh,
                              const bool ctr_phys_mesh,
                              const bool prm_mesh,
                              const bool sol_prm_mesh,
                              const bool knt_prm_mesh,
                              vtkMultiBlockDataSet *const mb) const;


  /**
   * @brief Checks if a given file can be read.
   *
   * It check that the file exists, is not corrupted, has the right
   * permissions, etc.
   *
   * If the file is binary, check that igatools serialization capabilities
   * are activated, elsewhere, if it is an ascii file, checks that
   * the XML igatools capabilities are activated.
   *
   * If the file can not be read, a @ref ExcVtkError exception is thrown,
   * that will be captured by @ref IgatoolsParaviewReader and shown
   * in the ParaView log window.
   */
  static void check_file(const std::string &file_name);

  /** @name Methods for querying and setting information of the VTK grids.*/
  ///@{

  /**
   * Returns the number of physical domains.
   * @return Number of physical domains.
   */
  Size get_number_physical_domains() const;

  /**
   * Returns the number of parametric domains.
   * @return Number of parametric domains.
   */
  Size get_number_parametric_domains() const;

  /**
   * @brief Returns the name of a physical domain identified by an
   * @p id number.
   *
   * Given the @p id number associated to the VTK grid, it returns its
   * physical domains' name.
   *
   * @param[in] Id of the VTK grid.
   * @return Name of the domain.
   */
  const char *get_physical_domain_name(const Index &id) const;

  /**
   * @brief Returns the name of a parametric domain identified by an
   * @p id number.
   *
   * Given the @p id number associated to the VTK grid, it returns its
   * parametric domains' name.
   *
   * @param[in] Id of the VTK grid.
   * @return Name of the domain.
   */
  const char *get_parametric_domain_name(const Index &id) const;

  /**
   * @brief Returns the status of a physical domain.
   *
   * Given the name of the physical domain, it returns the status (active
   * or not).
   *
   * @param[in] name Name of the domain.
   * @return true if the domain is active, false elsewhere.
   */
  bool get_physical_domain_status(const std::string &name) const;

  /**
   * @brief Returns the status of a parametric domain.
   *
   * Given the name of the parametric domain, it returns the status (active
   * or not).
   *
   * @param[in] name Name of the domain.
   * @return true if the domain is active, false elsewhere.
   */
  bool get_parametric_domain_status(const std::string &name) const;

  /**
   * @brief Sets the status of a physical domain.
   *
   * Given the name of the physical domain, it sets the status (active
   * or not).
   *
   * @param[in] name Name of the domain.
   * @param[in] status Status to be set in the domain.
   */
  void set_physical_domain_status(const std::string &name,
                                  const bool status);

  /**
   * @brief Sets the status of a parametric domain.
   *
   * Given the name of the parametric domain, it sets the status (active
   * or not).
   *
   * @param[in] name Name of the domain.
   * @param[in] status Status to be set in the domain.
   */
  void set_parametric_domain_status(const std::string &name,
                                    const bool status);

private:
  /**
   * @brief Retrieves the number of domains in a given grid container.
   *
   * @param[in] grid_container Container to be checked.
   * @return Number of domains.
   */
  static Size get_number_domains(const VtkIgaGridsContainer_ grid_container);

  /**
   * @brief Retrieves the number of \a active domains in a given grid container.
   *
   * @param[in] grid_container Container in which to look for the active
   * domains.
   * @return Number of domains.
   */
  static Size get_number_active_domains(const VtkIgaGridsContainer_ grid_container);

  /**
   * @brief Retrieves the name of a domain identified by its @p id
   * and contained in @p grid_container.
   *
   * @param[in] grid_container Container in which to look for the domain
   * with the given @p id.
   * @param[i] id Id number of the domain to be looked for.
   * @return Name of the domain.
   */
  static const char *get_domain_name(const VtkIgaGridsContainer_ grid_container,
                                     const Index &id);

  /**
   * @brief Retrieves the status of a domain identified by its @p name
   * and contained in @p grid_container.
   *
   * @param[in] grid_container Container in which to look for the domain
   * with the given @p name.
   * @param[i] name Name of the domain to be looked for.
   * @return Status (active or inactive) of the domain.
   */
  static bool get_domain_status(const VtkIgaGridsContainer_ grid_container,
                                const std::string &name);

  /**
   * @brief Sets the status of a domain identified by its @p name
   * and contained in @p grid_container.
   *
   * @param[in] grid_container Container in which to look for the domain
   * with the given @p name.
   * @param[i] name Name of the domain to be looked for.
   * @param[in] status Status (active or inactive) to be set.
   */
  static void set_domain_status(const VtkIgaGridsContainer_ grid_container,
                                const std::string &name,
                                const bool status);

  /**
   * @brief Creates the VTK grids of the solid meshes of the
   * domains and includes them into the multiblock data set.
   *
   * Iterates along @p grid_container and creates all the VTK grids
   * of the solid meshes of the domains that are active and inserts
   * them into the multiblock data set.
   *
   * The grids are created by calling to the @ref VtkIgaGrid methods,
   * that will call to the @ref VtkIgaSolidGrid helper class.
   *
   * @param[in] grid_container Container for creating the required grids.
   * @param[out] mb Multiblock data set for inserting the created VTK grids.
   */
  static void set_solid_domains(const VtkIgaGridsContainer_ grid_container,
                                vtkMultiBlockDataSet *const mb);

  /**
   * @brief Creates the VTK grids of the knot meshes of the
   * domains and includes them into the multiblock data set.
   *
   * Iterates along @p grid_container and creates all the VTK grids
   * of the knot meshes of the domains that are active and inserts
   * them into the multiblock data set.
   *
   * The grids are created by calling to the @ref VtkIgaGrid methods,
   * that will call to the @ref VtkIgaKnotGrid helper class.
   *
   * @param[in] grid_container Container for creating the required grids.
   * @param[out] mb Multiblock data set for inserting the created VTK grids.
   */
  static void set_knot_domains(const VtkIgaGridsContainer_ grid_container,
                               vtkMultiBlockDataSet *const mb);

  /**
   * @brief Creates the VTK grids of the control polygon meshes of the
   * domains and includes them into the multiblock data set.
   *
   * Iterates along @p grid_container and creates all the VTK grids
   * of the control polygon meshes of the domains that are active and inserts
   * them into the multiblock data set.
   *
   * The grids are created by calling to the @ref VtkIgaGrid methods,
   * that will call to the @ref VtkIgaControlGrid helper class.
   *
   * @param[in] grid_container Container for creating the required grids.
   * @param[out] mb Multiblock data set for inserting the created VTK grids.
   */
  static void set_control_domains(const VtkIgaGridsContainer_ grid_container,
                                  vtkMultiBlockDataSet *const mb);

  ///@}


  /**
   * @brief Checks if a given file is written in binary format or not.
   * @param[in] file_name Name of the file to be checked.
   * @return true if the file is written in binary format, false elsewhere.
   */
  static bool is_file_binary(const std::string &file_name);

  /**
   * @brief Parses an igatools objects container from a given @p file_name.
   *
   * Before parsing the object container, it is checked if the file
   * is suitable by calling the @ref check_file method.
   *
   * If it is identify that the file has a binary format (by calling
   * the method @ref is_file_binary), the objects container is parsed
   * by means of the igatools serialization techniques.
   *
   * Otherwise, if the file is ascii, the objects container it is
   * parsed by means of the helper class @ref ObjectsContainerXMLReader.
   * This class parses objects container written in an ascii XML human
   * readable format. See @ref ObjectsContainerXMLReader and
   * @subpage in_out for further details.
   *
   * In the case that problems occurs during the parsing of the container,
   * errors will be thrown that must be caught by @ref IgatoolsParaViewReader.
   *
   * @param[in] file_name Name of the file to be parsed.
   * @return Objects container parsed.
   */
  static ObjContPtr_ parse_objects_container(const std::string &file_name);

  /**
   * @brief Extracts all the objects with invalid dimensions in an objects
   * container and returns a description of each o them in a string.
   *
   * The plugin is unable to visualize objects (grids, grid functions,
   * domains and functions) that have dimension 0 or greater than 3.
   *
   * Given an objects container, this methods identifies all the
   * that have an invalid dimension and writes a description of each of
   * them in a string. The description similar to:
   * <tt>
   *    "Grid<0>, Name: grid_0, ObjectId: 12"
   * </tt>
   *
   * @param[in] obj_container Objects container.
   * @return Vector containing the description in a string of all the
   * invalid dimension objects.
   */
  static SafeSTLVector<std::string>
  get_invalid_dimension_objects(const ObjContPtr_ obj_container);

  /**
   * @brief Given an @ref ObjectsContainer, creates a new complete.
   *
   * By receiving the objects container parsed from the input file,
   * this method creates a new objects container will all the dependencies
   * include (see @ref ObjectsContainer for further details) and
   * stores it in @ref objects_container_.
   *
   * E.g. in an objects container where only a @ref Domain has been
   * inserted, its associated @ref Grid will be also inserted
   * (this can occur when the objects container has been serialized to
   * a file that is loaded by the plugin).
   *
   * All the objects inserted in the new container @ref objects_container_
   * are inserted as @p const.
   *
   * For those objects with invalid dimensions (extracted by calling
   * @ref get_invalid_dimension_objects) a warning will be thrown
   * and shown in the ParaView log window.
   *
   * @param[in] obj_container Objects container to be completed.
   */
  void fill_objects_container(const ObjContPtr_ objs_container);

  /**
   * @brief Sets the names of all the objects inside the container.
   *
   * The already existing names are preserved.
   *
   * For those object missing a name, a new one will be created.
   * The new name will be related in its relationship with other objects
   * present in the container.
   * E.g. if there exits a domain called 'MyDomain', but its grid does not
   * have a name, the grid will be renamed as 'Grid of the Domain "MyDomain"'.
   *
   * When no possible relationship is found, the object will be renamed
   * based on its unique object id number, e.g. 'Function Id=14'.
   *
   *
   * In the case of repetition, the name is slightly modified for noting
   * this repetition. E.g. if two domains are called 'MyDomain', one of
   * them will be called 'MyDomain (2)'.
   */
  void set_container_names();

  /**
   * @brief Builds and stores all the @ref VtkIgaGrid associated to each
   * IGA domain.
   *
   * All the created @ref VtkIgaGrids are stored in two containers:
   * @ref  phys_grids_ (for physical domains) and @ref parm_grids_
   * for (parametric domains).
   */
  void build_vtk_iga_grids();

};

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_GRID_CONTAINER_H_
