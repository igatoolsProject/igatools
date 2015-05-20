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

#include "vtkUnstructuredGridAlgorithm.h"

class vtkIgatoolsReader : public vtkUnstructuredGridAlgorithm
{
public:

  vtkTypeMacro(vtkIgatoolsReader, vtkUnstructuredGridAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkIgatoolsReader *New();

  // Description:
  // Specify file name of the .iga file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Set the number of visualization points per direction for each Bezier
  // element.
  vtkSetVector3Macro (NumVisualizationPoints, int);
  vtkGetVectorMacro  (NumVisualizationPoints, int, 3);

  /*
   * Set/Get Grid type
   *  - 0: solid grid.
   *  - 1: control points mesh grid.
   *  - 2: unit mesh grid.
   */
  vtkSetMacro(GridType, int);
  vtkGetMacro(GridType,int);

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
   * Check number of visualization points.
   */
  void check_number_visualization_points ();

  /*
   * Grid type.
   */
  int GridType;

  /*
   * File name variable.
   */
  char* FileName = NULL;
  int   NumVisualizationPoints[3];
};


#endif // VTK_IGATOOLS_READER_H_