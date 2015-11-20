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

#ifndef IGATOOLS_READER_H_
#define IGATOOLS_READER_H_

#include <vtkMultiBlockDataSetAlgorithm.h>

#include <igatools/base/config.h>
#include <igatools/utils/tensor_size.h>

#include <paraview_plugin/types.h>

/** Forward declarations. */

namespace iga
{
template <int> class TensorSize;
struct VtkGridInformation;
struct VtkControlGridInformation;

class FunctionsContainer;

class VtkIgaGridGeneratorContParm;
class VtkIgaGridGeneratorContPhys;
}

/**
 * @brief This is main class called from the ParaView xml plugin.
 *
 * Once a igatools file (containing the serialization of a \ref FunctionsContainer
 * class instance) has been selected, it receives some information from the
 * ParaView GUI (the file name itself, information for building the vtk
 * geometries, etc) and creates and returns the vtk grids packed into
 * MultiBlockDataSet.
 *
 * Each mapping, i.e. the function
 * \f$ \mathbf{F} \colon \hat{\Omega}\in\mathbb{R}^{\text{dim}} \to
 * \Omega\in\mathbb{R}^{\text{dim}+\text{codim}} \f$
 * representing the geometry, is represented in ParaView by a vtk grid.
 *
 * For each mapping, different functions
 * \f$ \mathbf{g} \colon \hat{\Omega}\in\mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{range}\times\text{rank}} \f$
 * can be associated to it, being represented in ParaView throw vtk point data.
 *
 * The plugin it is also capable of visualizing the parametric representation
 * of the mapping throw the use of \ref IdentityFunction
 * \f$ \mathbf{Id} \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}} \f$
 * representing the geometry, is represented in ParaView by a vtk grid.
 *
 * If the options set in the ParaView GUI are updated, only needed grids
 * are recomputed.
 *
 *
 * @author P. Antolin, 2015.
 *
 * @ingroup paraview_plugin
 */

class VTK_EXPORT IgatoolsParaViewReader : public vtkMultiBlockDataSetAlgorithm
{

private:

  /** @name Types definition. */
  ///@{

  /** Type of the class itself. */
  typedef IgatoolsParaViewReader Self_;

  /** Type for shared pointer of \ref vtkGridInformation.  */
  typedef std::shared_ptr<iga::VtkGridInformation> GridInfoPtr_;

  /** Type for shared pointer of \ref vtkControlGridInformation.  */
  typedef std::shared_ptr<iga::VtkControlGridInformation> ControlGridInfoPtr_;

  /** Type for shared pointer of \ref VtkGridGeneratorContPhys.  */
  typedef std::shared_ptr<iga::VtkIgaGridGeneratorContPhys> PhysGenPtr_;

  /** Type for shared pointer of \ref VtkGridGeneratorContParm.  */
  typedef std::shared_ptr<iga::VtkIgaGridGeneratorContParm> ParmGenPtr_;

  ///@}



  /** @name Constructor and destructor. */
  ///@{

  /** Default constructor. */
  IgatoolsParaViewReader();


  /** Destructor. Defined by default. */
  ~IgatoolsParaViewReader() = default;

  /** Copy constructor. Not allowed to be used. */
  IgatoolsParaViewReader(const Self_ &);

  /** Move constructor. Not allowed to be used. */
  IgatoolsParaViewReader(const Self_ &&);

  ///@}


  /** @name Assignment operator. */
  ///@{

  /** Copy assignment operator. Not allowed to be used. */
  Self_ operator=(const Self_ &) = delete;

  /** Move assignment operator. Not allowed to be used. */
  Self_ operator=(const Self_ &&) = delete;

  ///@}


protected:

  /** @name ParaView plugin requested methods. */
  ///@{


  /**
   * This method is called by the vtk superclass and is the one
   * in charge of building the vtk objects (in this case, a
   * vtkMultiBlockDataSet) and store it in @p outputVector.
   *
   * The returned vtkMultiBlockDataSet is a block of blocks.
   * That is, the returned block is composed by two blocks: one for the
   * physical mappings and another one for the parametric ones.
   *
   * The physical mappings block is itself a vtkMultiBlockDataSet composed
   * by the blocks of the solid, knot mesh and control mesh representation.
   * Every of these objects themselves are a vtkMultiBlockDataSet containing
   * one grid for every active mapping.
   *
   * In the same way, the parametric mappings blocks, is composed by the
   * solid and knot mesh blocks, that are also vtkMultiBlockDataSet.
   *
   * This method internally updates all the grid generators, that, if needed,
   * compute the new grids.
   */
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *) override final;


  /**
   * This vtk method in thought to be in charge of testing the input file
   * and obtaining minimal information from it (i.e. number of grids, types,
   * etc).
   *
   * However, due to the structure and the current capabilities of igatools
   * I/O functions, this method actually full parses the input file,
   * retrieves the func_container_, and creates the phys_gen_ and parm_gen_
   * generators.
   * Due to the fact that isogeometric geometries and results should not
   * be very memory consuming (compared to other methods), this actions should
   * be too expensive.
   *
   */
  virtual int RequestInformation(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) override final;

  virtual int FillOutputPortInformation(int port, vtkInformation *info) override final;

public:

  /** Required for ParaView plugin. */
  static IgatoolsParaViewReader *New();

  /** Required for ParaView plugin. */
  vtkTypeMacro(IgatoolsParaViewReader, vtkMultiBlockDataSetAlgorithm);

  /**
   * Prints the information of the class in the stream @p using the
   * indentation @p indent
   */
  void PrintSelf(ostream &os, vtkIndent indent) override;


  /**
   * Return the status (active/inactive) of the physical mapping defined by
   * its @p name
   */
  int GetPhysGeomArrayStatus(const char *name);

  /**
   * Sets the @p status (active/inactive) for the current ParaView
   * visualization of the physical mapping defined by its @p name
   */
  void SetPhysGeomArrayStatus(const char *name, int status);

  /**
   * Returns the number of physical mappings that can be visualized.
   */
  int GetNumberOfPhysGeomArrays();

  /**
   * Returns the name of the physical mapping with number @p index
   */
  const char *GetPhysGeomArrayName(int index);


  /**
   * Return the status (active/inactive) of the parametric mapping defined by
   * its @p name
   */
  int GetParmGeomArrayStatus(const char *name);

  /**
   * Sets the @p status (active/inactive) for the current ParaView
   * visualization of the parametric mapping defined by its @p name
   */
  void SetParmGeomArrayStatus(const char *name, int status);

  /**
   * Returns the number of parametric mappings that can be visualized.
   */
  int GetNumberOfParmGeomArrays();

  /**
   * Returns the name of the physical mapping with number @p index
   */
  const char *GetParmGeomArrayName(int index);


  /**
   * Sets the @p name of the igatools input file containing a serialization
   * of a \ref FunctionsContainer class instance.
   */
  virtual void SetFileName(const char *name);

  /**
   * Sets the number of vtk cells in each direction (i.e. @p arg1, @p arg2,
   * @p arg3) that will be used for visualizing each Bezier element
   * for the solid representation of a physical mapping.
   *
   * For dimension 2, the number of cells will be (@p arg1, @p arg2) and
   * for dim = 1, @p arg1 cells will be used in each direction.
   */
  virtual void SetNumVisualizationElementsPhysicalSolid(int arg1, int arg2, int arg3);


  /**
   * Sets the number of vtk cells in each direction (i.e. @p arg1, @p arg2,
   * @p arg3) that will be used for visualizing each Bezier element
   * for the solid representation of a parametric mapping.
   *
   * For dimension 2, the number of cells will be (@p arg1, @p arg2) and
   * for dim = 1, @p arg1 cells will be used in each direction.
   */
  virtual void SetNumVisualizationElementsParametricSolid(int arg1, int arg2, int arg3);

  /**
   * Sets the number of vtk cells in each direction (i.e. @p arg1, @p arg2,
   * @p arg3) that will be used for visualizing each Bezier element
   * for the knot mesh representation of a physical mapping.
   *
   * For dimension 2, the number of cells will be (@p arg1, @p arg2) and
   * for dim = 1, @p arg1 cells will be used in each direction.
   */
  virtual void SetNumVisualizationElementsPhysicalKnot(int arg1, int arg2, int arg3);

  /**
   * Sets the number of vtk cells in each direction (i.e. @p arg1, @p arg2,
   * @p arg3) that will be used for visualizing each Bezier element
   * for the knot mesh representation of a parametric mapping.
   *
   * For dimension 2, the number of cells will be (@p arg1, @p arg2) and
   * for dim = 1, @p arg1 cells will be used in each direction.
   */
  virtual void SetNumVisualizationElementsParametricKnot(int arg1, int arg2, int arg3);

  /*
   * Sets the \ref vtkGridType for the solid representation of the
   * physical mappings. The possible different values of @p arg are:
   *  - 0: unstructured vtk grid : quadratic cells.
   *  - 1: unstructured vtk grid : linear cells.
   *  - 2: structured vtk grid.
   */
  virtual void SetGridTypePhysicalSolid(int arg);

  /*
   * Sets the \ref vtkGridType for the solid representation of the
   * parametric mappings. The possible different values of @p arg are:
   *  - 0: unstructured vtk grid : quadratic cells.
   *  - 1: unstructured vtk grid : linear cells.
   *  - 2: structured vtk grid.
   */
  virtual void SetGridTypeParametricSolid(int arg);

  /*
   * Sets the \ref vtkGridType for the knot mesh representation of the
   * physical mappings. The possible different values of @p arg are:
   *  - 0: unstructured vtk grid : quadratic cells.
   *  - 1: unstructured vtk grid : linear cells.
   */
  virtual void SetGridTypePhysicalKnot(int arg);

  /*
   * Sets the \ref vtkGridType for the knot mesh representation of the
   * parametric mappings. The possible different values of @p arg are:
   *  - 0: unstructured vtk grid : quadratic cells.
   *  - 1: unstructured vtk grid : linear cells.
   */
  virtual void SetGridTypeParametricKnot(int arg);

  /*
   * Sets the \ref vtkGridType for the control mesh representation of the
   * physical mappings. The possible different values of @p arg are:
   *  - 1: unstructured vtk grid : linear cells.
   *  - 2: structured grid.
   */
  virtual void SetGridTypePhysicalControl(int arg);


  /*
   * Sets active/inactive the creation of the solid representation of the
   * physical mappings. The two possible values of @p arg are:
   *  - true:  create the vtk grids.
   *  - false: do not create vtk grids.
   */
  virtual void SetSolidMeshPhysical(bool arg);

  /*
   * Sets active/inactive the creation of the solid representation of the
   * parametric mappings. The two possible values of @p arg are:
   *  - true:  create the vtk grids.
   *  - false: do not create vtk grids.
   */
  virtual void SetSolidMeshParametric(bool arg);

  /*
   * Sets active/inactive the creation of the control mesh representation of
   * the physical mappings. The two possible values of @p arg are:
   *  - true:  create the vtk grids.
   *  - false: do not create vtk grids.
   */
  virtual void SetControlMeshPhysical(bool arg);

  /*
   * Sets active/inactive the creation of the knot mesh representation of the
   * physical mappings. The two possible values of @p arg are:
   *  - true:  create the vtk grids.
   *  - false: do not create vtk grids.
   */
  virtual void SetKnotMeshPhysical(bool arg);

  /*
   * Sets active/inactive the creation of the knot mesh representation of the
   * parametric mappings. The two possible values of @p arg are:
   *  - true:  create the vtk grids.
   *  - false: do not create vtk grids.
   */
  virtual void SetKnotMeshParametric(bool arg);

  /*
   * Sets active/inactive the creation of the of the physical mapping
   * representations. The two possible values of @p arg are:
   *  - true:  create the vtk grids.
   *  - false: do not create vtk grids.
   */
  virtual void SetPhysicalMesh(bool arg);

  /*
   * Sets active/inactive the creation of the of the parametric mapping
   * representations. The two possible values of @p arg are:
   *  - true:  create the vtk grids.
   *  - false: do not create vtk grids.
   */
  virtual void SetParametricMesh(bool arg);

  /*
   * Test whether the file with the given @p name exists and can be read by
   * this reader. If the file is Ok, returns 1, elsewhere, returns 0.
   */
  int CanReadFile(const char *name);

  ///@}


private:

  /*
   * This method is used by others (all the SetNumVisualizationElementsXXX)
   * to fill the passed by argument variable @p arr, with the values @p arg1,
   * @p arg2 and @ arg2.
   *
   * It also uses the @p name of the variable and the @p mesh_type for
   * debugging purposes.
   */
  void set_num_vis_elements(int arg1, int arg2, int arg3,
                            const char *const name,
                            const char *const mesh_type,
                            iga::TensorSize<3> &arr);

  /*
   * This method is used by others (all the SetGridTypeXXX)
   * to fill the passed by argument vtkGridType @p type, as a function
   * of the value @p arg1. The values of @p arg1 correspond to:
   *
   *  - 0: unstructured vtk grid : quadratic cells.
   *  - 1: unstructured vtk grid : linear cells.
   *  - 2: structured vtk grid.
   *
   * It also uses the @p name of the variable for debugging purposes.
   */
  void set_grid_type(int arg1,
                     const char *const name,
                     iga::vtkGridType &type);

  /** vtkGridType for the solid representation of the physical mappings. */
  iga::vtkGridType phys_sol_grid_type_ = iga::vtkGridType::None;

  /** vtkGridType for the solid representation of the parametric mappings. */
  iga::vtkGridType parm_sol_grid_type_ = iga::vtkGridType::None;

  /** vtkGridType for the knot mesh representation of the physical mappings. */
  iga::vtkGridType phys_knt_grid_type_ = iga::vtkGridType::None;

  /** vtkGridType for the knot mesh representation of the parametric mappings. */
  iga::vtkGridType parm_knt_grid_type_ = iga::vtkGridType::None;

  /** vtkGridType for the control mesh representation of the physical mappings. */
  iga::vtkGridType phys_ctr_grid_type_ = iga::vtkGridType::None;

  /**
   *  Flag for determining if the solid representation of the physical
   *  mapping grids must be created (true) or not (false).
   */
  bool create_sol_mesh_phys_ = false;

  /**
   *  Flag for determining if the solid representation of the parametric
   *  mapping grids must be created (true) or not (false).
   */
  bool create_sol_mesh_parm_ = false;

  /**
   *  Flag for determining if the knot mesh representation of the
   *  physical grids mappings must be created (true) or not (false).
   */
  bool create_knt_mesh_phys_ = false;

  /**
   *  Flag for determining if the knot mesh representation of the
   *  parametric mapping grids must be created (true) or not (false).
   */
  bool create_knt_mesh_parm_ = false;

  /**
   *  Flag for determining if the control mesh representation of the
   *  physical mapping grids must be created (true) or not (false).
   */
  bool create_ctr_mesh_phys_ = false;

  /**
   *  Flag for determining if the physical mapping grids must be
   *  created (true) or not (false).
   */
  bool create_physical_mesh_ = false;

  /**
   *  Flag for determining if the parametric mapping grids must be
   *  created (true) or not (false).
   */
  bool create_parametric_mesh_ = false;

  /**
   * Flag for indicating if the file must be parsed.
   *
   * @Note: the file must be parsed every time a new input file is set.
   */
  bool parse_file_ = true;

  /**
   * Full name (including path) for the input file.
   */
  char *file_name_ = NULL;

  /**
   * Sets the number of vtk cells in each direction that will be used for
   * visualizing each Bezier element for the solid representation of the
   * the physical mappings.
   *
   * For dimension 3, all the components of the tensor will be used.
   *
   * For dimension 2, the number of cells will be
   * (@p n_vis_elem_phys_solid_[0], @p @p n_vis_elem_phys_solid_[1]) and
   * for dim = 1, @p n_vis_elem_phys_solid_[0] cells will be used in each
   * direction.
   */
  iga::TensorSize<3> n_vis_elem_phys_solid_;

  /**
   * Sets the number of vtk cells in each direction that will be used for
   * visualizing each Bezier element for the solid representation of the
   * the parametric mappings.
   *
   * For dimension 3, all the components of the tensor will be used.
   *
   * For dimension 2, the number of cells will be
   * (@p n_vis_elem_phys_solid_[0], @p @p n_vis_elem_phys_solid_[1]) and
   * for dim = 1, @p n_vis_elem_phys_solid_[0] cells will be used in each
   * direction.
   */
  iga::TensorSize<3> n_vis_elem_parm_solid_;

  /**
   * Sets the number of vtk cells in each direction that will be used for
   * visualizing each Bezier element for the knot mesh representation of the
   * the physical mappings.
   *
   * For dimension 3, all the components of the tensor will be used.
   *
   * For dimension 2, the number of cells will be
   * (@p n_vis_elem_phys_solid_[0], @p @p n_vis_elem_phys_solid_[1]) and
   * for dim = 1, @p n_vis_elem_phys_solid_[0] cells will be used in each
   * direction.
   */
  iga::TensorSize<3> n_vis_elem_phys_knot_ ;

  /**
   * Sets the number of vtk cells in each direction that will be used for
   * visualizing each Bezier element for the knot mesh representation of the
   * the parametric mappings.
   *
   * For dimension 3, all the components of the tensor will be used.
   *
   * For dimension 2, the number of cells will be
   * (@p n_vis_elem_phys_solid_[0], @p @p n_vis_elem_phys_solid_[1]) and
   * for dim = 1, @p n_vis_elem_phys_solid_[0] cells will be used in each
   * direction.
   */
  iga::TensorSize<3> n_vis_elem_parm_knot_;


  /**
   * Igatools functions container used for storing the mapping functions
   * and their associated field data (@see FunctionsContainer).
   */
  std::shared_ptr<iga::FunctionsContainer> funcs_container_;


  /**
   * Generator for creating physical mappings (and associated field
   * data) representations throw a tree of vtkMultiBlockDataSet.
   */
  PhysGenPtr_ phys_gen_;

  /**
   * Generator for creating parametric mappings (and associated field
   * data) representations throw a tree of vtkMultiBlockDataSet.
   */
  ParmGenPtr_ parm_gen_;

  /**
   * Parses the input file.
   * It returns 1 if everything worked well, elsewhere returns 0.
   */
  int parse_file();

  /**
   * Updates the information related to the vtk grids for the grids generators.
   */
  void update_grid_info();

  /**
   * This is a temporary class for filling a functions container while
   * the serialization does not work.
   */
  template <int dim>
  void create_geometries();

  /**
   * Creates the tree of vtkMultiBlockDataSet of the output, and calls the
   * grid generators for filling the tree.
   */
  int create_grids(vtkMultiBlockDataSet *const mb);
};

#endif // IGATOOLS_READER_H_
