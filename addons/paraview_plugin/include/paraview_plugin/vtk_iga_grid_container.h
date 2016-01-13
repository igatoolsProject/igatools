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

/**
 * @brief To be defined.
 *
 * @todo To be defined.
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
template <class Domain> class VtkIgaGrid;
struct VtkGridInformation;
struct VtkControlGridInformation;

class VtkIgaGridContainer
{
private:

  /**
   * Self type.
   */
  typedef VtkIgaGridContainer Self_;

  /**
   * Self shared pointer type.
   */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /**
   * Grid generator type for @p dim and @p codim.
   */
  template <class Dmn>
  using GridGenPtr_ = std::shared_ptr<VtkIgaGrid<Dmn>>;

  typedef SafeSTLArray<int, 3> NumCells_;

  template <class T>
  struct IsInValidDim :
    boost::mpl::or_<
  boost::mpl::bool_<(T::dim < 1)>,
  boost::mpl::bool_<(T::dim > 3)>>
                                {};

  template <class T>
  struct IsInValidSpaceDim :
    boost::mpl::or_<
  boost::mpl::bool_<(T::space_dim < 1)>,
  boost::mpl::bool_<(T::space_dim > 3)>>
                                      {};

  template <class T>
  struct IsInValidDomain :
    boost::mpl::or_<
    IsInValidDim<T>,
    IsInValidSpaceDim<T>>
  {};

  template <class T>
  struct IsInValidFunction :
    boost::mpl::or_<
    IsInValidDim<T>,
    IsInValidSpaceDim<T>>
  {};

  template< class T >
  struct as_fusion_vector_shared_ptr
  {
    /**
     * This functor transform a <tt>boost::mpl::vector</tt> of types into a
     *  <tt>boost::fusion::vector</tt> of <tt>shared_ptr</tt>s of the types.
     */

    typedef typename boost::fusion::result_of::as_vector<
    typename boost::mpl::transform<T, std::shared_ptr<boost::mpl::_1>>::type>::type type;
  };

  /**
   * Valid domains.
   */
  using ValidDomains_ = boost::mpl::remove_if<
                        InstantiatedTypes::Domains,
                        boost::mpl::lambda< IsInValidDomain< boost::mpl::_1 > >::type
                        >::type;

  /**
   * Valid grids.
   */
  using ValidGrids_ = boost::mpl::remove_if<
                      InstantiatedTypes::Grids,
                      boost::mpl::lambda< IsInValidDim< boost::mpl::_1 > >::type
                      >::type;

  /**
   * Valid grids functions.
   */
  using ValidGridFuncs_ = typename boost::mpl::remove_if<
                          InstantiatedTypes::GridFunctions,
                          typename boost::mpl::lambda< IsInValidDim< boost::mpl::_1 > >::type
                          >::type;

  /**
   * Valid functions.
   */
  using ValidFunctions_ = typename boost::mpl::remove_if<
                          InstantiatedTypes::Functions,
                          typename boost::mpl::lambda< IsInValidDomain< boost::mpl::_1 > >::type
                          >::type;

  /**
   * Invalid domains.
   */
  using InvalidDomains_ = boost::mpl::remove_if<
                        InstantiatedTypes::Domains,
                        boost::mpl::lambda< boost::mpl::not_<IsInValidDomain< boost::mpl::_1 >> >::type
                        >::type;

  /**
   * Invalid grids.
   */
  using InvalidGrids_ = boost::mpl::remove_if<
                      InstantiatedTypes::Grids,
                      boost::mpl::lambda< boost::mpl::not_<IsInValidDim< boost::mpl::_1 >> >::type
                      >::type;

  /**
   * Invalid grids functions.
   */
  using InvalidGridFuncs_ = typename boost::mpl::remove_if<
                          InstantiatedTypes::GridFunctions,
                          boost::mpl::lambda< boost::mpl::not_<IsInValidDim< boost::mpl::_1 >> >::type
                          >::type;

  /**
   * Invalid functions.
   */
  using InvalidFunctions_ = typename boost::mpl::remove_if<
                          InstantiatedTypes::Functions,
                          boost::mpl::lambda< boost::mpl::not_<IsInValidDomain< boost::mpl::_1 >> >::type
                          >::type;

  using GridPtrs_ = as_fusion_vector_shared_ptr<ValidGrids_>::type;
  using GridFuncPtrs_ = as_fusion_vector_shared_ptr<ValidGridFuncs_>::type;
  using DomainPtrs_ = as_fusion_vector_shared_ptr<ValidDomains_>::type;
  using FunctionPtrs_ = as_fusion_vector_shared_ptr<ValidFunctions_>::type;

  using InvalidGridPtrs_ = as_fusion_vector_shared_ptr<InvalidGrids_>::type;
  using InvalidGridFuncPtrs_ = as_fusion_vector_shared_ptr<InvalidGridFuncs_>::type;
  using InvalidDomainPtrs_ = as_fusion_vector_shared_ptr<InvalidDomains_>::type;
  using InvalidFunctionPtrs_ = as_fusion_vector_shared_ptr<InvalidFunctions_>::type;

  template< class T >
  struct as_fusion_vector_const_shared_ptr
  {
  private:
    template <class S>
    using Pair_ = boost::fusion::pair<S, SafeSTLVector<GridGenPtr_<S>>>;

  public:
    typedef typename boost::fusion::result_of::as_map<
    typename boost::mpl::transform<T, Pair_<boost::mpl::_1>>::type>::type type;
  };



  template <class T, class Domain>
  struct IsInValidFunctionForDomain_ :
    boost::mpl::or_<
  boost::mpl::bool_<(T::dim != Domain::dim)>,
  boost::mpl::bool_<(T::space_dim != Domain::space_dim)>>
                                                       {};

  template <class T, class Domain>
  struct IsInValidGridFunctionForDomain_ :
    boost::mpl::or_<
  boost::mpl::bool_<(T::dim != Domain::dim)>,
  boost::mpl::bool_<(Domain::space_dim != Domain::dim)>>
                                                      {};

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
   * Basic container for grid generators.
   */
  using GridGensContainer_ = as_fusion_vector_const_shared_ptr<ValidDomains_>::type;

  /// Grid information shared pointer type.
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /// Objects container shared pointer type.
  typedef std::shared_ptr<ObjectsContainer> ObjContPtr_;


  /// Control grid information shared pointer type.
  typedef std::shared_ptr<VtkControlGridInformation> ControlGridInfoPtr_;

  /**
   * Void constructor.
   */
  VtkIgaGridContainer();

  /**
   * Constructor.
   */
  VtkIgaGridContainer(const ObjContPtr_ objs_container,
                      const GridInfoPtr_ phys_solid_info,
                      const GridInfoPtr_ phys_knot_info,
                      const ControlGridInfoPtr_ phys_control_info,
                      const GridInfoPtr_ parm_solid_info,
                      const GridInfoPtr_ parm_knot_info,
                      const std::string &file_name);

public:

  /**
   * TODO: to document.
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
   * TODO: to document.
   */
  static SelfPtr_ create_void();

  /**
   * TODO: to document.
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

public:
  /**
   * TODO: to document.
   */
  void create_multiblock_grid(const bool phys_mesh,
                              const bool sol_phys_mesh,
                              const bool knt_phys_mesh,
                              const bool ctr_phys_mesh,
                              const bool prm_mesh,
                              const bool sol_prm_mesh,
                              const bool knt_prm_mesh,
                              vtkMultiBlockDataSet *const mb) const;




  static void check_file (const std::string &file_name);

private:
  static bool is_file_binary(const std::string &file_name);
  static ObjContPtr_ parse_objects_container(const std::string &file_name);

  static SafeSTLVector<std::string>
  get_invalid_dimension_objects(const ObjContPtr_ obj_container);

  void fill_objects_container(const ObjContPtr_ objs_container);

  void set_names();

  void build_generators();

  template <class Domain>
  void insert_generator(const std::shared_ptr<const Domain> domain);

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

  /// Collection of physical grid generators.
  GridGensContainer_ phys_generators_;

  /// Collection of parametric grid generators.
  GridGensContainer_ parm_generators_;

  /// File name.
  std::string file_name_;

  /**
   * Container for numbering the generators included in the container.
   * Each entry of the vector is a tuple, whose components are:
   *   - 1st: global id of domain (global in the igatools framework).
   *   - 2nd: name associated to the domain function.
   *   - 3rd: flag for indicating if it is active or not.
   *   - 4rd: flag for indicating if it is an ig grid function or not.
   */
  SafeSTLVector<std::tuple<Index, std::string, bool, bool>> generators_numbering_;

public:

  /**
   * Returns the number of physical grids.
   */
  Size get_number_physical_grids() const;

  /**
   * Returns the number of parametric grids.
   */
  Size get_number_parametric_grids() const;

  /**
   * Returns the name of the @p id physical grid.
   * */
  const char *get_physical_grid_name(const Index &id) const;

  /**
   * Returns the name of the @p id parametric grid.
   */
  const char *get_parametric_grid_name(const Index &id) const;

  /**
   * Returns the status (active/inactive) of the physical grid named
   * @p name.
   */
  bool get_physical_grid_status(const std::string &name) const;

  /**
   * Returns the status (active/inactive) of the parametric grid named
   * @p name.
   */
  bool get_parametric_grid_status(const std::string &name) const;

  /**
   * Set the @p status (active/inactive) of the physical grid named
   * @p name.
   */
  void set_physical_grid_status(const std::string &name,
                                const bool status);

  /**
   * Set the @p status (active/inactive) of the parametric grid named
   * @p name.
   */
  void set_parametric_grid_status(const std::string &name,
                                  const bool status);

private:

  /**
   * TODO: to document.
   */
  static Size get_number_grids(const GridGensContainer_ generators);

  /**
   * TODO: to document.
   */
  static Size get_number_active_grids(const GridGensContainer_ generators);

  /**
   * TODO: to document.
   */
  static const char *get_grid_name(const GridGensContainer_ generators,
                                   const Index &id);

  /**
   * TODO: to document.
   */
  static bool get_grid_status(const GridGensContainer_ generators,
                              const std::string &name);

  /**
   * TODO: to document.
   */
  static void set_grid_status(const GridGensContainer_ generators,
                              const std::string &name,
                              const bool status);

  /**
   * TODO: to document.
   */
  static void set_solid_grids(const GridGensContainer_ generators,
                              vtkMultiBlockDataSet *const mb);

  /**
   * TODO: to document.
   */
  static void set_knot_grids(const GridGensContainer_ generators,
                             vtkMultiBlockDataSet *const mb);

  /**
   * TODO: to document.
   */
  static void set_control_grids(const GridGensContainer_ generators,
                                vtkMultiBlockDataSet *const mb);

};

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_GRID_CONTAINER_H_
