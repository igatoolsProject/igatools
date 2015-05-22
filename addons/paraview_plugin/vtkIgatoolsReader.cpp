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
  this->DebugOn ();
#endif

  this->NumVisualizationPoints[0] = 2;
  this->NumVisualizationPoints[1] = 2;
  this->NumVisualizationPoints[2] = 2;

  this->SetNumberOfInputPorts(0); // No vtk input, this is not a filter.
  this->SetNumberOfOutputPorts(1); // Just one output.
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

  this->check_number_visualization_points ();

  vtkInformation* info = outputVector->GetInformationObject(0);

  vtkDataObject* output = info->Get(vtkDataObject::DATA_OBJECT());
  vtkMultiBlockDataSet* mb =  vtkMultiBlockDataSet::SafeDownCast(output);

  if (!mb)
    return 0;

  this->SetErrorCode(vtkErrorCode::NoError);

  // TODO: to make it working here.
//   // Two blocks: a block for identity maps and another for physical maps.
//   mb->SetNumberOfBlocks(2);
//   mb->SetBlock(0, vtkSmartPointer<vtkMultiBlockDataSet>::New ());
//   mb->SetBlock(1, vtkSmartPointer<vtkMultiBlockDataSet>::New ());
// 
//   unsigned int index = 0;
//   mb->GetMetaData(index)->Set(vtkCompositeDataSet::NAME(), "Identity maps");
//   ++index;
//   mb->GetMetaData(index)->Set(vtkCompositeDataSet::NAME(), "Physical maps");

  return 1;


#if 0

 // This may be wrong if geometry is not cached. It is
 // update below.
 int numBlocks = static_cast<int>(this->Internal->Blocks.size());

 // Don't read the geometry if we already have it!
 if ( numBlocks == 0 )
   {
   FILE* xyzFp;
   if ( this->CheckGeometryFile(xyzFp) != VTK_OK)
     {
     return 0;
     }

   if ( this->ReadGeometryHeader(xyzFp) != VTK_OK )
     {
     vtkErrorMacro("Error reading geometry file.");
     fclose(xyzFp);
     return 0;
     }

   // Update from the value in the file.
   numBlocks = static_cast<int>(this->Internal->Blocks.size());

   for(i=0; i<numBlocks; i++)
     {

     // Read the geometry of this grid.
     this->SkipByteCount(xyzFp);

     vtkStructuredGrid* nthOutput = this->Internal->Blocks[i];
     int dims[3];
     nthOutput->GetDimensions(dims);
     vtkDataArray* pointArray = this->NewFloatArray();
     pointArray->SetNumberOfComponents(3);
     pointArray->SetNumberOfTuples( dims[0]*dims[1]*dims[2] );

     vtkPoints* points = vtkPoints::New();
     points->SetData(pointArray);
     pointArray->Delete();
     nthOutput->SetPoints(points);
     points->Delete();
     if ( this->ReadVector(xyzFp,
                           dims[0]*dims[1]*dims[2],
                           this->Internal->NumberOfDimensions,
                           pointArray) == 0)
       {
       vtkErrorMacro("Encountered premature end-of-file while reading "
                     "the geometry file (or the file is corrupt).");
       this->SetErrorCode(vtkErrorCode::PrematureEndOfFileError);
       fclose(xyzFp);
       return 0;
       }

     if (this->Internal->IBlanking)
       {
       int* ib = (int*)malloc(dims[0]*dims[1]*dims[2]*sizeof(int));
       if ( this->ReadIntBlock(xyzFp, dims[0]*dims[1]*dims[2], ib) == 0)
         {
         vtkErrorMacro("Encountered premature end-of-file while reading "
                       "the q file (or the file is corrupt).");
         this->SetErrorCode(vtkErrorCode::PrematureEndOfFileError);
         free(ib);
         fclose(xyzFp);
         return 0;
         }

       vtkIntArray* iblank = vtkIntArray::New();
       iblank->SetName("IBlank");
       iblank->SetVoidArray(ib, dims[0]*dims[1]*dims[2], 0);
       nthOutput->GetPointData()->AddArray(iblank);
       iblank->Delete();

       vtkUnsignedCharArray* visibility = vtkUnsignedCharArray::New();
       visibility->SetNumberOfComponents(1);
       visibility->SetNumberOfTuples( nthOutput->GetNumberOfCells() );
       visibility->SetName("Visibility");
       nthOutput->SetCellVisibilityArray(visibility);
       nthOutput->GetCellData()->AddArray(visibility);
       vtkIdList* ids = vtkIdList::New();
       ids->SetNumberOfIds(8);
       vtkIdType numCells = nthOutput->GetNumberOfCells();
       for (vtkIdType cellId=0; cellId<numCells; cellId++)
          {
          nthOutput->GetCellPoints(cellId, ids);
          vtkIdType numIds = ids->GetNumberOfIds();
          char visible = 1;
          for (vtkIdType ptIdx=0; ptIdx<numIds; ptIdx++)
            {
            if (ib[ids->GetId(ptIdx)] == 0)
              {
              visible = 0;
              break;
              }
            }
          visibility->SetValue(cellId, visible);
          }
        ids->Delete();
        }
      this->SkipByteCount(xyzFp);
      }

    fclose(xyzFp);
    }

  // Now read the solution.
  if (this->QFileName && this->QFileName[0] != '\0')
    {
    FILE* qFp;
    if ( this->CheckSolutionFile(qFp) != VTK_OK)
      {
      return 0;
      }

    int nq, nqc, isOverflow;
    if ( this->ReadQHeader(qFp, true, nq, nqc, isOverflow) != VTK_OK )
      {
      fclose(qFp);
      return 0;
      }

    for(i=0; i<numBlocks; i++)
      {
      vtkStructuredGrid* nthOutput = this->Internal->Blocks[i];


      // Save the properties first
      vtkDataArray* properties = this->NewFloatArray();
      properties->SetName("Properties");

      int numProperties = 4;
      int count = this->SkipByteCount(qFp);
      // We have a byte count to tell us how many Q values to
      // read. If this is more that 4, this is probably an Overflow
      // file.
      if (isOverflow)
        {
        // We take 4 bytes because there is an int there that
        // we will throw away
        numProperties = (count-4) / this->Internal->Precision + 1;
        }
      properties->SetNumberOfTuples(numProperties);

      // Read fsmach, alpha, re, time;
      if ( this->ReadScalar(qFp, 4, properties) == 0)
        {
        vtkErrorMacro("Encountered premature end-of-file while reading "
                      "the q file (or the file is corrupt).");
        this->SetErrorCode(vtkErrorCode::PrematureEndOfFileError);
        fclose(qFp);
        properties->Delete();
        return 0;
        }

      if (isOverflow)
        {
        // We create a dummy array to use with ReadScalar
        vtkDataArray* dummyArray = properties->NewInstance();
        dummyArray->SetVoidArray(properties->GetVoidPointer(4), 3, 1);

        // Read GAMINF, BETA, TINF
        if ( this->ReadScalar(qFp, 3, dummyArray) == 0)
          {
          vtkErrorMacro("Encountered premature end-of-file while reading "
                        "the q file (or the file is corrupt).");
          this->SetErrorCode(vtkErrorCode::PrematureEndOfFileError);
          fclose(qFp);
          properties->Delete();
          return 0;
          }

        // igam is an int
        int igam;
        this->ReadIntBlock(qFp, 1, &igam);
        properties->SetTuple1(7, igam);

        dummyArray->SetVoidArray(properties->GetVoidPointer(8), 3, 1);
        // Read the rest of properties
        if ( this->ReadScalar(qFp, numProperties - 8, dummyArray) == 0)
          {
          vtkErrorMacro("Encountered premature end-of-file while reading "
                        "the q file (or the file is corrupt).");
          this->SetErrorCode(vtkErrorCode::PrematureEndOfFileError);
          fclose(qFp);
          properties->Delete();
          return 0;
          }
        dummyArray->Delete();
        }

      nthOutput->GetFieldData()->AddArray(properties);
      properties->Delete();
      this->SkipByteCount(qFp);

      int dims[3];
      nthOutput->GetDimensions(dims);

      this->SkipByteCount(qFp);

      vtkDataArray* density = this->NewFloatArray();
      density->SetNumberOfComponents(1);
      density->SetNumberOfTuples( dims[0]*dims[1]*dims[2] );
      density->SetName("Density");
      if ( this->ReadScalar(qFp, dims[0]*dims[1]*dims[2], density) == 0)
        {
        vtkErrorMacro("Encountered premature end-of-file while reading "
                      "the q file (or the file is corrupt).");
        this->SetErrorCode(vtkErrorCode::PrematureEndOfFileError);
        fclose(qFp);
        density->Delete();
        return 0;
        }
      nthOutput->GetPointData()->AddArray(density);
      density->Delete();

      vtkDataArray* momentum = this->NewFloatArray();
      momentum->SetNumberOfComponents(3);
      momentum->SetNumberOfTuples( dims[0]*dims[1]*dims[2] );
      momentum->SetName("Momentum");
      if ( this->ReadVector(qFp,
                            dims[0]*dims[1]*dims[2],
                            this->Internal->NumberOfDimensions,
                            momentum) == 0)
        {
        vtkErrorMacro("Encountered premature end-of-file while reading "
                      "the q file (or the file is corrupt).");
        this->SetErrorCode(vtkErrorCode::PrematureEndOfFileError);
        fclose(qFp);
        momentum->Delete();
        return 0;
        }
      nthOutput->GetPointData()->AddArray(momentum);
      momentum->Delete();

      vtkDataArray* se = this->NewFloatArray();
      se->SetNumberOfComponents(1);
      se->SetNumberOfTuples( dims[0]*dims[1]*dims[2] );
      se->SetName("StagnationEnergy");
      if (this->ReadScalar(qFp, dims[0]*dims[1]*dims[2], se) == 0)
        {
        vtkErrorMacro("Encountered premature end-of-file while reading "
                      "the q file (or the file is corrupt).");
        fclose(qFp);
        se->Delete();
        return 0;
        }
      nthOutput->GetPointData()->AddArray(se);
      se->Delete();

      if (isOverflow)
        {
        if(nq >= 6) 
          {
          vtkDataArray* gamma = this->NewFloatArray();
          gamma->SetNumberOfComponents(1);
          gamma->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
          gamma->SetName("Gamma");
          if (this->ReadScalar(qFp, dims[0]*dims[1]*dims[2], gamma) == 0)
            {
            vtkErrorMacro("Encountered premature end-of-file while reading "
                          "the q file (or the file is corrupt).");
            fclose(qFp);
            gamma->Delete();
            return 0;
            }
          nthOutput->GetPointData()->AddArray(gamma);
          gamma->Delete();
        } 

        char res[100];
        // Read species and turbulence variables for overflow q files
        for(int j=0; j<nqc; j++)
          {
          vtkDataArray* temp = this->NewFloatArray();
          temp->SetNumberOfComponents(1);
          temp->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
          int k = j+1;
          sprintf(res, "Species Density #%d", k);
          temp->SetName(res);
          if (this->ReadScalar(qFp, dims[0]*dims[1]*dims[2], temp) == 0)
            {
            vtkErrorMacro("Encountered premature end-of-file while reading "
                          "the q file (or the file is corrupt).");
            fclose(qFp);
            temp->Delete();
            return 0;
            }
          nthOutput->GetPointData()->AddArray(temp);
          temp->Delete();
          }
        float d, r;
        for(int v=0; v<nqc; v++)
          {
          vtkDataArray* rat = this->NewFloatArray();
          sprintf(res, "Species Density #%d", v+1);
          vtkPointData* outputPD = nthOutput->GetPointData();
          vtkDataArray* spec = outputPD->GetArray(res);
          vtkDataArray* dens = outputPD->GetArray("Density");
          rat->SetNumberOfComponents(1);
          rat->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
          sprintf(res, "Spec Dens #%d / rho", v+1);
          rat->SetName(res);
          for(int w=0; w<dims[0]*dims[1]*dims[2]; w++)
            {
            r = dens->GetComponent(w,0);
            r = (r != 0.0 ? r : 1.0);
            d = spec->GetComponent(w,0);
            rat->SetTuple1(w, d/r);
            }
          nthOutput->GetPointData()->AddArray(rat);
          rat->Delete();
          }
        for(int a=0; a<nq-6-nqc; a++)
          {
          vtkDataArray* temp = this->NewFloatArray();
          temp->SetNumberOfComponents(1);
          temp->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
          int k = a+1;
          sprintf(res, "Turb Field Quant #%d", k);
          temp->SetName(res);
          if (this->ReadScalar(qFp, dims[0]*dims[1]*dims[2], temp) == 0)
            {
            vtkErrorMacro("Encountered premature end-of-file while reading "
                          "the q file (or the file is corrupt).");
            fclose(qFp);
            temp->Delete();
            return 0;
            }
          nthOutput->GetPointData()->AddArray(temp);
          temp->Delete();
          }
        }

      this->SkipByteCount(qFp);

      if ( this->FunctionList->GetNumberOfTuples() > 0 )
        {
        int fnum;
        for (int tup=0; tup < this->FunctionList->GetNumberOfTuples(); tup++)
          {
          if ( (fnum=this->FunctionList->GetValue(tup)) >= 0 )
            {
            this->MapFunction(fnum, nthOutput);
            }
          }
        }
      this->AssignAttribute(this->ScalarFunctionNumber, nthOutput,
                            vtkDataSetAttributes::SCALARS);
      this->AssignAttribute(this->VectorFunctionNumber, nthOutput,
                            vtkDataSetAttributes::VECTORS);
      }
    fclose(qFp);
    }

  // Now read the functions.
  if (this->FunctionFileName && this->FunctionFileName[0] != '\0')
    {
    FILE* fFp;
    if ( this->CheckFunctionFile(fFp) != VTK_OK)
      {
      return 0;
      }

    std::vector<int> nFunctions(numBlocks);
    if ( this->ReadFunctionHeader(fFp, &nFunctions[0]) != VTK_OK )
      {
      fclose(fFp);
      return 0;
      }

    for(i=0; i<numBlocks; i++)
      {
      vtkStructuredGrid* nthOutput = this->Internal->Blocks[i];
      int dims[3];
      nthOutput->GetDimensions(dims);

      this->SkipByteCount(fFp);

      for (int j=0; j<nFunctions[i]; j++)
        {
        vtkDataArray* functionArray = this->NewFloatArray();
        functionArray->SetNumberOfTuples( dims[0]*dims[1]*dims[2] );
        char functionName[20];
        sprintf(functionName, "Function%d", j);
        functionArray->SetName(functionName);
        if (this->ReadScalar(fFp, dims[0]*dims[1]*dims[2], functionArray) == 0)
          {
          vtkErrorMacro("Encountered premature end-of-file while reading "
                        "the function file (or the file is corrupt).");
          fclose(fFp);
          functionArray->Delete();
          return 0;
          }
        nthOutput->GetPointData()->AddArray(functionArray);
        functionArray->Delete();
        }

      this->SkipByteCount(fFp);
      }
    }
#endif

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

  iga_vtk_.set_file (file_name, file_path);
  iga_vtk_.set_number_visualization_points (this->GetNumVisualizationPoints ());

  vtkInformation* info = outputVector->GetInformationObject(0);
  vtkDataObject* output = info->Get(vtkDataObject::DATA_OBJECT());
  vtkMultiBlockDataSet* mb =  vtkMultiBlockDataSet::SafeDownCast(output);

  // Two blocks: a block for identity maps and another for physical maps.
  mb->SetNumberOfBlocks(2);
  mb->SetBlock(0, vtkSmartPointer<vtkMultiBlockDataSet>::New ());
  mb->SetBlock(1, vtkSmartPointer<vtkMultiBlockDataSet>::New ());

  unsigned int index = 0;
  mb->GetMetaData(index)->Set(vtkCompositeDataSet::NAME(), "Identity maps");
  ++index;
  mb->GetMetaData(index)->Set(vtkCompositeDataSet::NAME(), "Physical maps");

  iga_vtk_.parse_file ();
  iga_vtk_.generate_vtk_grids (mb);
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
