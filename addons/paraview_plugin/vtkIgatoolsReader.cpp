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

#include <vtkIgatoolsReader.h>

#include <paraview_plugin/iga_grid_generator.h>

#include <vtkObjectFactory.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>

#include <vtkSMSourceProxy.h>

#include <sys/stat.h>

using std::string;


vtkStandardNewMacro(vtkIgatoolsReader);

vtkIgatoolsReader::vtkIgatoolsReader()
{
  this->DebugOn ();
  vtkDebugMacro (<< "vtkIgatoolsReader constructor begin\n");
  this->NumVisualizationPoints[0] = 2;
  this->NumVisualizationPoints[1] = 2;
  this->NumVisualizationPoints[2] = 2;
  this->GridType = 0;
  this->SetNumberOfInputPorts(0); // No vtk input, this is not a filter.
  this->SetNumberOfOutputPorts(1); // Just one output.
  vtkDebugMacro (<< "vtkIgatoolsReader constructor end\n");
//   vtkSMSourceProxy* poSourceProxy = vtkSMSourceProxy::SafeDownCast(this->poSourceProxy);
//   poSourceProxy->InvokeCommand("MyFunction");
}


int vtkIgatoolsReader::RequestInformation(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  vtkDebugMacro (<< "vtkIgatoolsReader RequestInformation begin\n");
  vtkDebugMacro (<< "vtkIgatoolsReader RequestInformation: file name"
                 << FileName << "\n");
  vtkDebugMacro (<< "vtkIgatoolsReader RequestInformation: number of "
                 << "visualization points: "
                 << NumVisualizationPoints[0] << " "
                 << NumVisualizationPoints[1] << " "
                 << NumVisualizationPoints[2] << "\n");
  vtkDebugMacro (<< "vtkIgatoolsReader RequestInformation end\n");

  if (this->GetGridType() != 1) // Solid grid or unit mesh.
  {
    this->check_number_visualization_points ();
    return 1;
  }
  else // Control mesh grid.
    return 1;


};

void
vtkIgatoolsReader::
check_number_visualization_points ()
{
  if (NumVisualizationPoints[0] < 2 || NumVisualizationPoints[1] < 2 ||
      NumVisualizationPoints[2] < 2)
  {
    int* num_points = this->GetNumVisualizationPoints ();
    vtkWarningMacro (<< "In vtkIgatoolsReader invalid specified number of visualization points "
                     << "per Bezier element (" 
                     << NumVisualizationPoints[0] << ", "
                     << NumVisualizationPoints[1] << ", "
                     << NumVisualizationPoints[2] << "). "
                     << "All the values must be >=2.\n"
                     << "The number of points was automatically set to ("
                     << (num_points[0] > 1 ? num_points[0] : 2) << ", "
                     << (num_points[1] > 1 ? num_points[1] : 2) << ", "
                     << (num_points[2] > 1 ? num_points[2] : 2) << ").\n");
    if (num_points[0] < 2) num_points[0] = 2;
    if (num_points[1] < 2) num_points[1] = 2;
    if (num_points[2] < 2) num_points[2] = 2;

    this->SetNumVisualizationPoints (num_points);
  }
};



int vtkIgatoolsReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  string file_name;
  string file_path;
  this->get_file_and_path (file_name, file_path);
  const auto grid_generator = IGAVTKGridGenerator::create (file_name, file_path, this->GetNumVisualizationPoints ());

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  if (this->GetGridType() == 0) // Solid grid.
    return grid_generator->fill_solid_output (outInfo);
  else if (this->GetGridType() == 2) // Unit mesh
    return grid_generator->fill_identity_output (outInfo);
  else // Control mesh grid.
    return grid_generator->fill_control_mesh_output (outInfo);
}

void vtkIgatoolsReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: "
      << (this->FileName ? this->FileName : "(none)") << "\n";
};



void
vtkIgatoolsReader::
get_file_and_path (string& file_name, string& file_path)
{
  file_name = string (this->GetFileName());
  file_path = file_name.substr (file_name.find_last_of ("\\/") + 1,
                                file_name.size ());
};



int 
vtkIgatoolsReader::
CanReadFile(const char *name)
{
  // First make sure the file exists.  This prevents an empty file
  // from being created on older compilers.
  struct stat fs;
  if(stat(name, &fs) != 0)
  {
    return 0;
  }

  // TODO: To check here the file extension and some additonal info of the final
  // for checking its suitableness.

  return 1;
};
