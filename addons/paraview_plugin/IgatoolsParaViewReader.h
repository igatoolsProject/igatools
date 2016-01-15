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

#ifndef __IGATOOLS_READER_H_
#define __IGATOOLS_READER_H_

#include <vtkMultiBlockDataSetAlgorithm.h>

#include <paraview_plugin/vtk_iga_types.h>
#include <igatools/utils/safe_stl_array.h>

class vtkObjectBase;

namespace iga
{
    template <class T, int d> class SafeSTLArray;
    namespace paraview_plugin
    {
        class VtkIgaGridContainer;
    }
}

/**
 * @brief Main class of the igatools-ParaView plugin.
 *
 * The main purpose of this class is receiving and sending information
 * to ParaView. Thus, it receives the input (the input file and VTK grid
 * options) and information for activating or de-activating each
 * IGA domain, and returns the VTK multiblock grid itself, and information
 * about the number of domains, their names, etc.
 *
 * The plugin receives from the ParaView GUI an input file name
 * and group of options for building the information for building the VTK
 * geometries. With this information creates and returns the VTK grids packed
 * into MultiBlockDataSet.
 *
 * The interaction between this class and ParaView is performed basically
 * throw two methods:
 *  - @ref RequestInformation: once it is called, the plugin creates new
 *    @ref paraview_plugin::VtkIgaGridContainer objects
 *    (@ref iga_grid_container_), that will parse the input file,
 *    and create a @ref paraview_plugin::VtkIgaGrid associated to each
 *    domain present in the input file.
 *  - @ref RequestData: when it is called, @ref iga_grid_container_
 *    is requested to create the VTK multi block object containing all
 *    the active domains.
 *
 * @see paraview_plugin::VtkIgaGrid
 * @see paraview_plugin::VtkIgaGridContainer
 * @see paraview_plugin::VtkGridInformation
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
class IgatoolsParaViewReader : public vtkMultiBlockDataSetAlgorithm
{

private:

  /// Type of the class itself.
  typedef IgatoolsParaViewReader Self_;

  /// Type for shared pointer of @ref VtkIgaGridContainer.
  typedef std::shared_ptr<iga::paraview_plugin::VtkIgaGridContainer> IgaGridContPtr_;

  /// Type for containing the number of cells per side of each Bezier element.
  typedef iga::SafeSTLArray<int, 3> NumCells_;


  /** @name Constructor and destructor. */
  ///@{
  /**
   * Constructor.
   */
  IgatoolsParaViewReader();


  /**
   * Default destructor.
   */
  ~IgatoolsParaViewReader() = default;

  /**
   * Copy constructor.
   * @warning Not allowed to be used.
   */
  IgatoolsParaViewReader(const Self_ &) = delete;

  /**
   * Move constructor.
   * @warning Not allowed to be used.
   */
  IgatoolsParaViewReader(const Self_ &&) = delete;

  ///@}


  /** @name Assignment operator. */
  ///@{

  /**
   * Copy assignment operator.
   * @warning Not allowed to be used.
   */
  Self_ operator=(const Self_ &) = delete;

  /**
   * Move assignment operator.
   * @warning Not allowed to be used.
   */
  Self_ operator=(const Self_ &&) = delete;

  ///@}


protected:

  /** @name ParaView plugin requested methods. */
  ///@{

  /**
   * @brief Creates the IGA geometries.
   *
   * This method is called by the VTK superclass and is the one
   * in charge of building the VTK objects (in this case, a
   * @ref vtkMultiBlockDataSet) and store it in @p output_vec.
   *
   * The returned @ref vtkMultiBlockDataSet is a block of blocks
   * (see @ref VtkIgaGrid for the details of its internal structure).
   *
   * Before creating the VTK grids, this method internally updates all
   * VTK grid information and communicates if to the class @ref VtkIgaGrid.
   *
   * @note Required by the ParaView plugin.
   *
   * @param[in] request Requested information.
   * @param[in] input_vec Passed input information
   * @param[out] output_vec Object receiving the output.
   * In this case receives nothing.
   * @return 1 if the process succeeded, 0 otherwise.
   */
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *) override final;


  /**
   * @brief Loads the plugin input file.
   *
   * This VTK method in thought to be in charge of testing the input file
   * and obtaining minimal information from it (i.e. number of grids, types,
   * etc).
   *
   * However, due to the structure and the current capabilities of igatools
   * I/O functions, this method actually full parses the input file,
   * retrieves an objects container and creates the @ref iga_grid_container_.
   *
   * Due to the fact that isogeometric geometries and results should not
   * be very memory consuming (compared to other methods), this actions
   * should not be too expensive.
   *
   * @note Required by the ParaView plugin.
   *
   * @param[in] request Requested information.
   * @param[in] input_vec Passed input information
   * @param[out] output_vec Object receiving the output.
   * In this case receives nothing.
   * @return 1 if the process succeeded, 0 otherwise.
   *
   */
  virtual int RequestInformation(vtkInformation *request,
                                 vtkInformationVector **input_vec,
                                 vtkInformationVector *output_vec) override final;

public:

  /**
   * @brief Retrieves the class name.
   * @note Required by the ParaView plugin.
   * @return Class name.
   */
  static IgatoolsParaViewReader *New();

private:
  /**
   * @brief Retrieves the class name.
   * @note Implementation extracted from <tt>common/vtkGetSet.h</tt>
   * @note Required by the ParaView plugin.
   * @return Class name.
   */
  virtual const char *GetClassNameInternal() const override final;

public:
  /**
   * @brief Returns true if the class is of the given type.
   * @note Implementation extracted from <tt>common/vtkGetSet.h</tt>
   * @note Required by the ParaView plugin.
   * @return true, if the current class is of the given type, false,
   * otherwise.
   */
  static int IsTypeOf(const char *type);

  /**
   * @brief Returns true if the current class is of the given type.
   * @note Implementation extracted from <tt>common/vtkGetSet.h</tt>
   * @note Required by the ParaView plugin.
   * @return true, if the current class is of the given type, false,
   * otherwise.
   */
  virtual int IsA(const char *type) override final;

  /**
   * @brief Perform a safe down cast of a given @ref vtkObjectBase
   * to the current type.
   *
   * @note Implementation extracted from <tt>common/vtkGetSet.h</tt>
   * @note Required by the ParaView plugin.
   *
   * @return Casted object.
   */
  static IgatoolsParaViewReader* SafeDownCast(vtkObjectBase *o);


  /**
   * @brief Creates a new instance of the class.
   *
   * @note Implementation extracted from <tt>common/vtkGetSet.h</tt>
   * @note Required by the ParaView plugin.
   *
   * @return New instance.
   */
  IgatoolsParaViewReader *NewInstance() const;

protected:

  /**
   * @brief Creates a new instance of the class down casted to a
   * @ref vtkObjectBase.
   *
   * @note Implementation extracted from <tt>common/vtkGetSet.h</tt>
   * @note Required by the ParaView plugin.
   *
   * @return New instance.
   */
  virtual vtkObjectBase *NewInstanceInternal() const override final;

public:

  /**
   * @brief Prints the information of the class in the stream @p using the
   * indentation @p indent
   * @param[in] os Stream for writing the information.
   * @param[in] indent Indentation to be used for the writing.
   */
  void PrintSelf(ostream &os, vtkIndent indent) override;


  /**
   * @brief Return the status (active/inactive) of the physical domain
   * defined by its @p name
   * @param[in] name Name of the domain to be checked.
   * @return Status of the domain (1 active, 0 inactive).
   */
  int GetPhysGeomArrayStatus(const char *name);

  /**
   * @brief Sets the @p status (active/inactive) for the current ParaView
   * visualization of the physical domain defined by its @p name
   * @param[in] name Name of the domain.
   * @param[in] status Status to be set.
   */
  void SetPhysGeomArrayStatus(const char *name, int status);

  /**
   * @brief Returns the number of physical domains.
   * @return Number of the physical domains.
   */
  int GetNumberOfPhysGeomArrays();

  /**
   * @brief Returns the name of the physical domain with number @p index
   *
   * @param[in] Number of the physical domain.
   */
  const char *GetPhysGeomArrayName(int index);


  /**
   * @brief Return the status (active/inactive) of the parametric domain
   * defined by its @p name
   * @param[in] name Name of the domain to be checked.
   * @return Status of the domain (1 active, 0 inactive).
   */
  int GetParmGeomArrayStatus(const char *name);

  /**
   * @brief Sets the @p status (active/inactive) for the current ParaView
   * visualization of the parametric domain defined by its @p name
   * @param[in] name Name of the domain.
   * @param[in] status Status to be set.
   */
  void SetParmGeomArrayStatus(const char *name, int status);

  /**
   * @brief Returns the number of parametric domains.
   * @return Number of the parametric domains.
   */
  int GetNumberOfParmGeomArrays();

  /**
   * @brief Returns the name of the physical domain with number @p index
   *
   * @param[in] Number of the parametric domain.
   */
  const char *GetParmGeomArrayName(int index);


  /**
   * @brief Sets the @p name of the igatools input file containing the
   * @ref ObjectsContainer class instance.
   *
   * @param[in] name Name of the input file.
   */
  virtual void SetFileName(const char *name);

  /**
   * @brief Sets the number of VTK cells in each direction (i.e. @p arg1, @p arg2,
   * @p arg3) that will be used for visualizing each Bezier element
   * for the solid representation of a physical domain.
   *
   * For dimension 2, the number of cells will be (@p arg1, @p arg2) and
   * for dim = 1, @p arg1 cells will be used in each direction.
   */
  virtual void SetNumVisualizationElementsPhysicalSolid(int arg1, int arg2, int arg3);


  /**
   * @brief Sets the number of VTK cells in each direction (i.e. @p arg1, @p arg2,
   * @p arg3) that will be used for visualizing each Bezier element
   * for the solid representation of a parametric domain.
   *
   * For dimension 2, the number of cells will be (@p arg1, @p arg2) and
   * for dim = 1, @p arg1 cells will be used in each direction.
   */
  virtual void SetNumVisualizationElementsParametricSolid(int arg1, int arg2, int arg3);

  /**
   * @brief Sets the number of VTK cells in each direction (i.e. @p arg1,
   * @p arg2, @p arg3) that will be used for visualizing each Bezier
   * element for the knot mesh representation of a physical domain.
   *
   * For dimension 2, the number of cells will be (@p arg1, @p arg2) and
   * for dim = 1, @p arg1 cells will be used in each direction.
   *
   * @param[in] arg1 Number of VTK cells in the first parametric
   * direction for each Bezier element.
   * @param[in] arg2 Number of VTK cells in the second parametric
   * direction for each Bezier element.
   * @param[in] arg3 Number of VTK cells in the third parametric
   * direction for each Bezier element.
   */
  virtual void SetNumVisualizationElementsPhysicalKnot(int arg1, int arg2, int arg3);

  /**
   * @brief Sets the number of VTK cells in each direction (i.e. @p arg1,
   * @p arg2, @p arg3) that will be used for visualizing each Bezier
   * element for the knot mesh representation of a parametric domain.
   *
   * For dimension 2, the number of cells will be (@p arg1, @p arg2) and
   * for dim = 1, @p arg1 cells will be used in each direction.
   *
   * @param[in] arg1 Number of VTK cells in the first parametric
   * direction for each Bezier element.
   * @param[in] arg2 Number of VTK cells in the second parametric
   * direction for each Bezier element.
   * @param[in] arg3 Number of VTK cells in the third parametric
   * direction for each Bezier element.
   */
  virtual void SetNumVisualizationElementsParametricKnot(int arg1, int arg2, int arg3);

  /**
   * @brief Sets the @ref VtkGridType for the solid mesh
   * representation of the physical domains.
   *
   * The possible different values of @p arg are:
   *  - 0: unstructured VTK grid : quadratic cells.
   *  - 1: unstructured VTK grid : linear cells.
   *  - 2: structured VTK grid.
   *
   *  @param[in] arg Value to be set.
   */
  virtual void SetGridTypePhysicalSolid(int arg);

  /**
   * @brief Sets the @ref VtkGridType for the solid mesh
   * representation of the parametric domains.
   *
   * The possible different values of @p arg are:
   *  - 0: unstructured VTK grid : quadratic cells.
   *  - 1: unstructured VTK grid : linear cells.
   *  - 2: structured VTK grid.
   *
   *  @param[in] arg Value to be set.
   */
  virtual void SetGridTypeParametricSolid(int arg);

  /**
   * @brief Sets the @ref VtkGridType for the knot mesh
   * representation of the physical domains.
   *
   * The possible different values of @p arg are:
   *  - 0: unstructured VTK grid : quadratic cells.
   *  - 1: unstructured VTK grid : linear cells.
   *
   *  @param[in] arg Value to be set.
   */
  virtual void SetGridTypePhysicalKnot(int arg);

  /**
   * @brief Sets the @ref VtkGridType for the knot mesh
   * representation of the parametric domains.
   *
   * The possible different values of @p arg are:
   *  - 0: unstructured VTK grid : quadratic cells.
   *  - 1: unstructured VTK grid : linear cells.
   *
   *  @param[in] arg Value to be set.
   */
  virtual void SetGridTypeParametricKnot(int arg);

  /**
   * @brief Sets the @ref VtkGridType for the control polygon mesh
   * representation of the physical domains.
   *
   * The possible different values of @p arg are:
   *  - 1: unstructured VTK grid : linear cells.
   *  - 2: structured grid.
   *
   *  @param[in] arg Value to be set.
   */
  virtual void SetGridTypePhysicalControl(int arg);


  /**
   * @brief Sets active/inactive the creation of the solid mesh
   * representation of the physical domains.
   *
   * The two possible values of @p arg are:
   *  - true:  create the VTK grids.
   *  - false: do not create VTK grids.
   *
   *  @param[in] arg Value to be set, true or false.
   */
  virtual void SetSolidMeshPhysical(bool arg);

  /**
   * @brief Sets active/inactive the creation of the solid mesh
   * representation of the parametric domains.
   *
   * The two possible values of @p arg are:
   *  - true:  create the VTK grids.
   *  - false: do not create VTK grids.
   *
   *  @param[in] arg Value to be set, true or false.
   */
  virtual void SetSolidMeshParametric(bool arg);

  /**
   * @brief Sets active/inactive the creation of the control polygon mesh
   * representation of the physical domains.
   *
   * The two possible values of @p arg are:
   *  - true:  create the VTK grids.
   *  - false: do not create VTK grids.
   *
   *  @param[in] arg Value to be set, true or false.
   */
  virtual void SetControlMeshPhysical(bool arg);

  /**
   * @brief Sets active/inactive the creation of the knot mesh
   * representation of the physical domains.
   *
   * The two possible values of @p arg are:
   *  - true:  create the VTK grids.
   *  - false: do not create VTK grids.
   *
   *  @param[in] arg Value to be set, true or false.
   */
  virtual void SetKnotMeshPhysical(bool arg);

  /**
   * @brief Sets active/inactive the creation of the knot mesh
   * representation of the parametric domains.
   *
   * The two possible values of @p arg are:
   *  - true:  create the VTK grids.
   *  - false: do not create VTK grids.
   *
   *  @param[in] arg Value to be set, true or false.
   */
  virtual void SetKnotMeshParametric(bool arg);

  /**
   * @brief Sets active/inactive the creation of the of the physical domain
   * representations.
   *
   * The two possible values of @p arg are:
   *  - true:  create the VTK grids.
   *  - false: do not create VTK grids.
   *
   *  @param[in] arg Value to be set, true or false.
   */
  virtual void SetPhysicalMesh(bool arg);

  /**
   * @brief Sets active/inactive the creation of the of the parametric domain
   * representations.
   *
   * The two possible values of @p arg are:
   *  - true:  create the VTK grids.
   *  - false: do not create VTK grids.
   *
   *  @param[in] arg Value to be set, true or false.
   */
  virtual void SetParametricMesh(bool arg);

  /**
   * Test whether the file with the given @p name exists and can be read by
   * this reader. If the file is Ok, returns 1, elsewhere, returns 0
   * and shows an error in the ParaView log window.
   *
   * @param[in] name Name of the file.
   */
  int CanReadFile(const char *name);

  ///@}


private:

  /**
   * @brief Sets the number of cells for each Bezier element to the
   * given @p arr.
   *
   * This method is used by others (all the @p SetNumVisualizationElementsXXX)
   * to fill the passed by argument variable @p arr, with the values @p arg1,
   * @p arg2 and @ arg3.
   *
   * @param[in] arg1 Number of VTK cells in the first parametric
   * direction for each Bezier element.
   * @param[in] arg2 Number of VTK cells in the second parametric
   * direction for each Bezier element.
   * @param[in] arg3 Number of VTK cells in the third parametric
   * direction for each Bezier element.
   * @param[in] name Name of the variable for being set, for debugging
   * purposes.
   * @param[in] mesh_type Type of the mesh, for debugging purposes.
   * @param[out] arr Variable where the number VTK cells are set.
   */
  void set_num_vis_elements(int arg1, int arg2, int arg3,
                            const char *const name,
                            const char *const mesh_type,
                            NumCells_ &arr);

  /**
   * @brief Sets the grid type to the given @p type.
   *
   * This method is used by others (all the @p SetGridTypeXXX)
   * to fill the passed by argument @ref VtkGridType @p type, as a function
   * of the value @p arg1. The values of @p arg1 correspond to:
   *
   *  - 0: unstructured VTK grid : quadratic cells.
   *  - 1: unstructured VTK grid : linear cells.
   *  - 2: structured VTK grid.
   *
   * @param[in] arg1 Value of the grid type to be set.
   * @param[in] name Name of the variable for being set, for debugging
   * purposes.
   * @param[out] type Variable where the value is set.
   */
  void set_grid_type(int arg1,
                     const char *const name,
                     iga::paraview_plugin::VtkGridType &type);

  /// VtkGridType for the solid representation of the physical domains.
  iga::paraview_plugin::VtkGridType phys_sol_grid_type_ = iga::paraview_plugin::VtkGridType::None;

  /// VtkGridType for the solid representation of the parametric domains.
  iga::paraview_plugin::VtkGridType parm_sol_grid_type_ = iga::paraview_plugin::VtkGridType::None;

  /// VtkGridType for the knot mesh representation of the physical domains.
  iga::paraview_plugin::VtkGridType phys_knt_grid_type_ = iga::paraview_plugin::VtkGridType::None;

  /// VtkGridType for the knot mesh representation of the parametric domains.
  iga::paraview_plugin::VtkGridType parm_knt_grid_type_ = iga::paraview_plugin::VtkGridType::None;

  /// VtkGridType for the control polygon mesh representation of the physical domains.
  iga::paraview_plugin::VtkGridType phys_ctr_grid_type_ = iga::paraview_plugin::VtkGridType::None;

  /**
   *  Flag for determining if the solid representation of the physical
   *  domains must be created (true) or not (false).
   */
  bool create_sol_mesh_phys_ = false;

  /**
   *  Flag for determining if the solid representation of the parametric
   *  domains must be created (true) or not (false).
   */
  bool create_sol_mesh_parm_ = false;

  /**
   *  Flag for determining if the knot mesh representation of the
   *  physical domains must be created (true) or not (false).
   */
  bool create_knt_mesh_phys_ = false;

  /**
   *  Flag for determining if the knot mesh representation of the
   *  parametric domains must be created (true) or not (false).
   */
  bool create_knt_mesh_parm_ = false;

  /**
   *  Flag for determining if the control polygon mesh representation of the
   *  physical domains must be created (true) or not (false).
   */
  bool create_ctr_mesh_phys_ = false;

  /**
   *  Flag for determining if the physical domains must be
   *  created (true) or not (false).
   */
  bool create_physical_mesh_ = false;

  /**
   *  Flag for determining if the parametric domains must be
   *  created (true) or not (false).
   */
  bool create_parametric_mesh_ = false;

  /**
   * Flag for indicating if the input file must be parsed.
   *
   * @Note: The input file must be parsed every time a new input file is set.
   */
  bool parse_input_file_ = true;

  /// Full name (including path) for the input file.
  char *file_name_ = NULL;

  /**
   * Number of VTK cells in each parametric direction for every Bezier
   * element for the solid mesh of the physical domain.
   */
  NumCells_ n_vis_elem_phys_solid_;

  /**
   * Number of VTK cells in each parametric direction for every Bezier
   * element for the solid mesh of the parametric domain.
   */
  NumCells_ n_vis_elem_parm_solid_;

  /**
   * Number of VTK cells in each parametric direction for every Bezier
   * element for the knot mesh of the physical domain.
   */
  NumCells_ n_vis_elem_phys_knot_;

  /**
   * Number of VTK cells in each parametric direction for every Bezier
   * element for the knot mesh of the parametric domain.
   */
  NumCells_ n_vis_elem_parm_knot_;

  /// Iga grid container for creating VTK grids.
  IgaGridContPtr_ iga_grid_container_;
};

#endif // __IGATOOLS_READER_H_
