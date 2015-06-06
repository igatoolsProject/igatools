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

#ifndef VTK_IGATOOLS_READER_H_
#define VTK_IGATOOLS_READER_H_

#include "vtkMultiBlockDataSetAlgorithm.h"

#include <paraview_plugin/iga_vtk.h>

class vtkStructuredGrid;

// class vtkIgatoolsReader : public vtkUnstructuredGridAlgorithm
class vtkIgatoolsReader : public vtkMultiBlockDataSetAlgorithm
{
public:

  vtkTypeMacro(vtkIgatoolsReader, vtkMultiBlockDataSetAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkIgatoolsReader *New();

  // Description:
  // Specify file name of the .iga file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Set the number of VTK visualization elements per direction for each Bezier
  // element.
  vtkSetVector3Macro (NumVisualizationElements, int);
  vtkGetVectorMacro  (NumVisualizationElements, int, 3);

  // Description:
  // Set the on/off the use of VTK quadratic elements for unstructured grids.
  vtkSetMacro (QuadraticElements, bool);
  vtkGetMacro (QuadraticElements, bool);

  /*
   * Set/Get Grid type
   *  - 0: unstructured grid.
   *  - 1: structured grid.
   */
  vtkSetMacro(GridType, int);
  vtkGetMacro(GridType, int);

  /*
   * Set/Get Control Mesh creation flag.
   *  - true:  create the mesh.
   *  - false: do not create the mesh.
   */
  vtkSetMacro(ControlMesh, bool);
  vtkGetMacro(ControlMesh, bool);

  /*
   * Set/Get Knot Mesh creation flag.
   *  - true:  create the mesh.
   *  - false: do not create the mesh.
   */
  vtkSetMacro(KnotMesh, bool);
  vtkGetMacro(KnotMesh, bool);

  /*
   * Set/Get Parametric Mesh creation flag.
   *  - true:  create the mesh.
   *  - false: do not create the mesh.
   */
  vtkSetMacro(ParametricMesh, bool);
  vtkGetMacro(ParametricMesh, bool);

  /*
   * Set/Get Physical Mesh creation flag.
   *  - true:  create the mesh.
   *  - false: do not create the mesh.
   */
  vtkSetMacro(PhysicalMesh, bool);
  vtkGetMacro(PhysicalMesh, bool);

protected:
  /*
   * Deleter
   */
  vtkIgatoolsReader();

  /*
   * Destructor.
   */
  ~vtkIgatoolsReader(){}

  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *) override final;

  virtual int RequestInformation(vtkInformation* request,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector) override final;

public:
  /*
   * Test whether the file with the given name exists and can be read by this
   * reader.
   */
  int CanReadFile(const char* name);

private:
  vtkIgatoolsReader(const vtkIgatoolsReader&) = delete;
  vtkIgatoolsReader(const vtkIgatoolsReader&&) = delete;
  void operator=(const vtkIgatoolsReader&) = delete;
  void operator=(const vtkIgatoolsReader&&) = delete;


  /*
   * Retrieves the file name and the path of the file.
   */
  void get_file_and_path (std::string& file_name, std::string& file_path);

  /*
   * Check number of visualization element.
   */
  void check_number_visualization_elements ();

  /*
   * Grid type.
   */
  int GridType;

  /*
   * Control mesh creation flag.
   */
  bool ControlMesh;

  /*
   * Knot mesh creation flag.
   */
  bool KnotMesh;

  /*
   * Parametric mesh creation flag.
   */
  bool ParametricMesh;

  /*
   * Physical mesh creation flag.
   */
  bool PhysicalMesh;

  /*
   * Quadratic elements flag.
   */
  bool QuadraticElements;

  /*
   * File name variable.
   */
  char* FileName = NULL;

  /*
   * Number of VTK visualization elements per direction per Bezier element.
   */
  int   NumVisualizationElements[3];

  /*
   * Iga vtk object.
   */
  iga::IGAVTK iga_vtk_;
};


#endif // VTK_IGATOOLS_READER_H_