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

#include <vtkObjectFactory.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSMSourceProxy.h>
#include <vtkErrorCode.h>
#include <vtkSmartPointer.h>


#include <sys/stat.h>

using std::string;


vtkStandardNewMacro(vtkIgatoolsReader);

vtkIgatoolsReader::vtkIgatoolsReader()
{
#ifndef NDEBUG
    this->DebugOn();
#endif

    this->NumVisualizationPoints[0] = 2;
    this->NumVisualizationPoints[1] = 2;
    this->NumVisualizationPoints[2] = 2;

    this->GridType = 0;
    this->ControlMesh = false;
    this->ParametricMesh = false;
    this->PhysicalMesh = false;

    this->SetNumberOfInputPorts(0); // No vtk input, this is not a filter.
    this->SetNumberOfOutputPorts(1); // Just one output.
}



int vtkIgatoolsReader::RequestInformation(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector)
{
    vtkDebugMacro(<< "vtkIgatoolsReader RequestInformation begin\n");
    vtkDebugMacro(<< "vtkIgatoolsReader RequestInformation: file name"
                  << FileName << "\n");
    vtkDebugMacro(<< "vtkIgatoolsReader RequestInformation: number of "
                  << "visualization points: "
                  << NumVisualizationPoints[0] << " "
                  << NumVisualizationPoints[1] << " "
                  << NumVisualizationPoints[2] << "\n");
    vtkDebugMacro(<< "vtkIgatoolsReader RequestInformation end\n");

    this->check_number_visualization_points();

    vtkInformation *info = outputVector->GetInformationObject(0);

    vtkDataObject *output = info->Get(vtkDataObject::DATA_OBJECT());
    vtkMultiBlockDataSet *mb =  vtkMultiBlockDataSet::SafeDownCast(output);

    if (!mb)
        return 0;

    this->SetErrorCode(vtkErrorCode::NoError);

    return 1;
};

void
vtkIgatoolsReader::
check_number_visualization_points()
{
    if (NumVisualizationPoints[0] < 2 || NumVisualizationPoints[1] < 2 ||
        NumVisualizationPoints[2] < 2)
    {
        int *num_points = this->GetNumVisualizationPoints();
        vtkWarningMacro(<< "In vtkIgatoolsReader invalid specified number of visualization points "
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

        this->SetNumVisualizationPoints(num_points);
    }
};



int vtkIgatoolsReader::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector)
{
    string file_name;
    string file_path;
    this->get_file_and_path(file_name, file_path);

    iga_vtk_.set_file(file_name, file_path);
    iga_vtk_.set_number_visualization_points(this->GetNumVisualizationPoints());

    vtkInformation *info = outputVector->GetInformationObject(0);
    vtkDataObject *output = info->Get(vtkDataObject::DATA_OBJECT());
    vtkMultiBlockDataSet *mb =  vtkMultiBlockDataSet::SafeDownCast(output);

    // Setting the blocks.
    int num_blocks = 0;
    if (this->GetControlMesh())
        ++num_blocks;
    if (this->GetParametricMesh())
        ++num_blocks;
    if (this->GetPhysicalMesh())
        ++num_blocks;

    mb->SetNumberOfBlocks(num_blocks);

    for (int i = 0; i < num_blocks; ++i)
        mb->SetBlock(i, vtkSmartPointer<vtkMultiBlockDataSet>::New());

    unsigned int index = 0;
    if (this->GetControlMesh())
        mb->GetMetaData(index++)->Set(vtkCompositeDataSet::NAME(), "Control mesh");
    if (this->GetParametricMesh())
        mb->GetMetaData(index++)->Set(vtkCompositeDataSet::NAME(), "Parametric mesh");
    if (this->GetPhysicalMesh())
        mb->GetMetaData(index++)->Set(vtkCompositeDataSet::NAME(), "Physical mesh");

    iga_vtk_.parse_file();
    iga_vtk_.generate_vtk_grids(this->GetGridType(),
                                this->GetControlMesh(),
                                this->GetParametricMesh(),
                                this->GetPhysicalMesh(),
                                mb);

    return 1;
}

void vtkIgatoolsReader::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);

    os << indent << "File Name: "
       << (this->FileName ? this->FileName : "(none)") << "\n";
};



void
vtkIgatoolsReader::
get_file_and_path(string &file_name, string &file_path)
{
    file_name = string(this->GetFileName());
    file_path = file_name.substr(file_name.find_last_of("\\/") + 1,
                                 file_name.size());
};



int
vtkIgatoolsReader::
CanReadFile(const char *name)
{
    // First make sure the file exists.  This prevents an empty file
    // from being created on older compilers.
    struct stat fs;
    if (stat(name, &fs) != 0)
    {
        return 0;
    }

    // TODO: To check here the file extension and some additonal info of the final
    // for checking its suitableness.

    return 1;
};

// 239
// 240 int vtkMultiBlockPLOT3DReader::CheckFile(FILE*& fp, const char* fname)
// 241 {
//   242   if (this->BinaryFile)
//   243     {
//     244     fp = fopen(fname, "rb");
//     245     }
//     246   else
//     247     {
//       248     fp = fopen(fname, "r");
//       249     }
//       250   if ( fp == NULL)
//       251     {
//         252     this->SetErrorCode(vtkErrorCode::FileNotFoundError);
//         253     vtkErrorMacro(<< "File: " << fname << " not found.");
//         254     return VTK_ERROR;
//         255     }
//         256   return VTK_OK;
//         257 }
//         258
//         259 int vtkMultiBlockPLOT3DReader::CheckGeometryFile(FILE*& xyzFp)
//         260 {
//           261   if ( this->XYZFileName == NULL || this->XYZFileName[0] == '\0'  )
//           262     {
//             263     this->SetErrorCode(vtkErrorCode::NoFileNameError);
//             264     vtkErrorMacro(<< "Must specify geometry file");
//             265     return VTK_ERROR;
//             266     }
//             267   return this->CheckFile(xyzFp, this->XYZFileName);
//             268 }
//             269
//             270 int vtkMultiBlockPLOT3DReader::CheckSolutionFile(FILE*& qFp)
//             271 {
//               272   if ( this->QFileName == NULL || this->QFileName[0] == '\0' )
//               273     {
//                 274     this->SetErrorCode(vtkErrorCode::NoFileNameError);
//                 275     vtkErrorMacro(<< "Must specify geometry file");
//                 276     return VTK_ERROR;
//                 277     }
//                 278   return this->CheckFile(qFp, this->QFileName);
//                 279 }
//                 280
//                 281 int vtkMultiBlockPLOT3DReader::CheckFunctionFile(FILE*& fFp)
//                 282 {
//                   283   if ( this->FunctionFileName == NULL || this->FunctionFileName[0] == '\0' )
//                   284     {
//                     285     this->SetErrorCode(vtkErrorCode::NoFileNameError);
//                     286     vtkErrorMacro(<< "Must specify geometry file");
//                     287     return VTK_ERROR;
//                     288     }
//                     289   return this->CheckFile(fFp, this->FunctionFileName);
//                     290 }
