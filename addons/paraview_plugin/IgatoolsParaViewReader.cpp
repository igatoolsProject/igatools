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

#include <IgatoolsParaViewReader.h>

#include <vtkObjectFactory.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>

#include <igatools/functions/identity_function.h>
#include <igatools/functions/functions_container.h>
#include <paraview_plugin/grid_generator_container.h>
#include <paraview_plugin/grid_information.h>

#include <sys/stat.h>


using std::get;
using std::string;
using std::shared_ptr;
using namespace iga;

vtkStandardNewMacro(IgatoolsParaViewReader);

IgatoolsParaViewReader::IgatoolsParaViewReader()
  :
  n_vis_elem_phys_solid_(1),
  n_vis_elem_parm_solid_(1),
  n_vis_elem_phys_knot_(1),
  n_vis_elem_parm_knot_(1),
  phys_gen_(PhysGenPtr_()),
  parm_gen_(ParmGenPtr_())
{
#ifndef NDEBUG
  this->DebugOn();
#else
  this->DebugOff();
#endif

  this->SetNumberOfInputPorts(0); // No vtk input, this is not a filter.
  this->SetNumberOfOutputPorts(1); // Just one output.

}



int IgatoolsParaViewReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *info = outputVector->GetInformationObject(0);

  vtkDataObject *output = info->Get(vtkDataObject::DATA_OBJECT());
  vtkMultiBlockDataSet *mb =  vtkMultiBlockDataSet::SafeDownCast(output);

  if (!mb)
    return 0;

  if (parse_file_)
    return this->parse_file();


  return 1;
}



int IgatoolsParaViewReader::FillOutputPortInformation(int port, vtkInformation *info)
{

  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}





int IgatoolsParaViewReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{


  vtkInformation *info = outputVector->GetInformationObject(0);
  vtkDataObject *output = info->Get(vtkDataObject::DATA_OBJECT());
  vtkMultiBlockDataSet *mb =  vtkMultiBlockDataSet::SafeDownCast(output);

  this->SetProgressText("Generating igatools geometries.");

  this->UpdateProgress(0.0);

  this->update_grid_info();
  return this->create_grids(mb);

}



void IgatoolsParaViewReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: "
     << (this->file_name_ ? this->file_name_ : "(none)") << "\n";
}



int
IgatoolsParaViewReader::
CanReadFile(const char *name)
{
  // TODO: it should be also determined here it the file is suitable,
  // and what kind of file it is.

  // Checking if the file exists.
  errno = 0;
  struct stat buffer;
  if (stat(name, &buffer) == -1) // ==0 ok; ==-1 error
  {
    string error_msg = string("Parsing file path ") + string(name) + string(" : ");
    if (errno == ENOENT)       // errno declared by include file errno.h
      error_msg += string("Path file does not exist, or path is an empty string.");
    else if (errno == ENOTDIR)
      error_msg += string("A component of the path is not a directory.");
    else if (errno == ELOOP)
      error_msg += string("Too many symbolic links encountered while traversing the path.");
    else if (errno == EACCES)
      error_msg += string("Permission denied.");
    else if (errno == ENAMETOOLONG)
      error_msg += string("File can not be read");
    else
      error_msg += string("An unknown problem was encountered.");
    vtkErrorMacro(<< error_msg);

    return 0;
  }
  else
    return 1;
}



int
IgatoolsParaViewReader::
parse_file()
{
  if (this->CanReadFile(file_name_) == 0)
    return 0;

  const auto file_name_str = string(file_name_);

  phys_gen_.reset();
  parm_gen_.reset();
  funcs_container_.reset();

//  Assert(false,ExcNotImplemented());

  try
  {

    this->SetProgressText("Parsing igatools file.");

    // TODO: this is a temporary solution for avoiding runtime error in deserialization.
    //   ifstream xml_istream(file_name_);
    //   IArchive xml_in(xml_istream);
    //   xml_in >> BOOST_SERIALIZATION_NVP (funcs_container_);
    //   xml_istream.close();
    //   Assert here if funcs_container_ is void.

    funcs_container_ = std::make_shared <FunctionsContainer> ();

    this->template create_geometries<1>();
    this->template create_geometries<2>();
    this->template create_geometries<3>();
    //*/

    // Physical solid grid.
    const auto phys_sol = VtkGridInformation::create
                          (n_vis_elem_phys_solid_, phys_sol_grid_type_);

    // Physical knot grid.
    const auto phys_knt = VtkGridInformation::create
                          (n_vis_elem_phys_knot_, phys_knt_grid_type_);

    // Physical control grid.
    const auto phys_ctr = VtkControlGridInformation::create
                          (phys_ctr_grid_type_ == vtkGridType::Structured);

    phys_gen_ = VtkIgaGridGeneratorContPhys::create
                (funcs_container_, phys_sol, phys_knt, phys_ctr);


    // Parametric solid grid.
    const auto parm_sol = VtkGridInformation::create
                          (n_vis_elem_parm_solid_, parm_sol_grid_type_);


    // Parametric knot grid.
    const auto parm_knt = VtkGridInformation::create
                          (n_vis_elem_parm_knot_, parm_knt_grid_type_);

    parm_gen_ = VtkIgaGridGeneratorContParm::create
                (funcs_container_, parm_sol, parm_knt);

    parse_file_ = false;
  }
  catch (std::exception &e)
  {
    vtkErrorMacro(<< e.what());

    phys_gen_.reset();
    parm_gen_.reset();
    funcs_container_.reset();

    return 0;
  }
  catch (...)
  {
    vtkErrorMacro(<< "An exception occurred when parsing file "
                  << "\"" << file_name_str << "\"");

    phys_gen_.reset();
    parm_gen_.reset();
    funcs_container_.reset();

    return 0;
  }

  return 1;
}



void
IgatoolsParaViewReader::
update_grid_info()
{
  // Physical solid grid.
  const auto phys_sol = VtkGridInformation::create
                        (n_vis_elem_phys_solid_, phys_sol_grid_type_);

  // Physical knot grid.
  const auto phys_knt = VtkGridInformation::create
                        (n_vis_elem_phys_knot_, phys_knt_grid_type_);

  // Physical control grid.
  const auto phys_ctr = VtkControlGridInformation::create
                        (phys_ctr_grid_type_ == vtkGridType::Structured);

  phys_gen_->update_physical(phys_sol, phys_knt, phys_ctr);


  // Parametric solid grid.
  const auto parm_sol = VtkGridInformation::create
                        (n_vis_elem_parm_solid_, parm_sol_grid_type_);


  // Parametric knot grid.
  const auto parm_knt = VtkGridInformation::create
                        (n_vis_elem_parm_knot_, parm_knt_grid_type_);

  parm_gen_->update_parametric(parm_sol, parm_knt);
}



int
IgatoolsParaViewReader::
create_grids(vtkMultiBlockDataSet *const mb)
{
  unsigned int num_blocks = create_physical_mesh_ + create_parametric_mesh_;

  if (num_blocks == 0)
  {
    vtkWarningMacro(<< "Neither physical nor parametric geometries are "
                    "active. No output produced.");

    return 1;
  }

  const unsigned int num_phys_blocks =
    create_sol_mesh_phys_ + create_knt_mesh_phys_ + create_ctr_mesh_phys_;

  const unsigned int num_parm_blocks =
    create_sol_mesh_parm_ + create_knt_mesh_parm_;

  const auto num_active_phys = phys_gen_->get_number_active_grids();
  const auto num_active_parm = parm_gen_->get_number_active_grids();

  bool new_create_physical_mesh = create_physical_mesh_;
  if (create_physical_mesh_ && (num_phys_blocks == 0 || num_active_phys == 0))
  {
    if (num_phys_blocks == 0)
    {
      vtkWarningMacro(<< "Physical geometries set active, but none grid type "
                      "(solid, knot, control) has been selected");
    }
    else
    {
      vtkWarningMacro(<< "Physical geometries set active, but no "
                      "geometry set active from the list.");
    }

    --num_blocks;

    new_create_physical_mesh = false;
  }


  bool new_create_parametric_mesh = create_parametric_mesh_;
  if (create_parametric_mesh_ && (num_parm_blocks == 0 || num_active_parm == 0))
  {
    if (num_parm_blocks == 0)
    {
      vtkWarningMacro(<< "Parametric geometries set active, but none grid type "
                      "(solid, knot) has been selected");
    }
    else
    {
      vtkWarningMacro(<< "Parametric geometries set active, but no "
                      "geometry set active from the list.");
    }

    --num_blocks;

    new_create_parametric_mesh = false;
  }

  if (num_blocks == 0)

    return 1;


  mb->SetNumberOfBlocks(num_blocks);

  // Creating blocks for the physical and parametric geometries.
  for (unsigned int i = 0; i < num_blocks; ++i)
    mb->SetBlock(i, vtkSmartPointer<vtkMultiBlockDataSet>::New());


  Size total_number_blocks = 0;
  if (new_create_physical_mesh)
  {
    if (create_sol_mesh_phys_)
      ++total_number_blocks;
    if (create_knt_mesh_phys_)
      ++total_number_blocks;
    if (create_ctr_mesh_phys_)
      ++total_number_blocks;
  }
  if (new_create_parametric_mesh)
  {
    if (create_sol_mesh_parm_)
      ++total_number_blocks;
    if (create_knt_mesh_parm_)
      ++total_number_blocks;
  }


  Index progress_index = 0;


  unsigned int block_index = 0;
  if (new_create_physical_mesh)
  {
    mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), "Physical mesh");


    vtkMultiBlockDataSet *const phys_block =
      vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index));

    phys_block->SetNumberOfBlocks(num_phys_blocks);


    Index subblock_index = 0;

    if (create_sol_mesh_phys_)
    {
      const auto solid_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      phys_block->SetBlock(subblock_index, solid_block);
      phys_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Solid mesh");
      phys_gen_->set_solid_grids(solid_block);

      this->UpdateProgress((++progress_index) / double (total_number_blocks));

      ++subblock_index;
    }

    if (create_knt_mesh_phys_)
    {
      const auto knot_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      phys_block->SetBlock(subblock_index, knot_block);
      phys_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Knot mesh");
      phys_gen_->set_knot_grids(knot_block);

      this->UpdateProgress((++progress_index) / double (total_number_blocks));

      ++subblock_index;
    }

    if (create_ctr_mesh_phys_)
    {
      const auto control_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      phys_block->SetBlock(subblock_index, control_block);
      phys_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Control mesh");
      phys_gen_->set_control_grids(control_block);

      this->UpdateProgress((++progress_index) / double (total_number_blocks));
    }


    ++block_index;
  } // create_physical_mesh


  if (new_create_parametric_mesh)
  {
    mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), "Parametric mesh");



    vtkMultiBlockDataSet *const parm_block =
      vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index));

    parm_block->SetNumberOfBlocks(num_parm_blocks);


    Index subblock_index = 0;

    if (create_sol_mesh_parm_)
    {
      const auto solid_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      parm_block->SetBlock(subblock_index, solid_block);
      parm_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Solid mesh");
      parm_gen_->set_solid_grids(solid_block);

      this->UpdateProgress((++progress_index) / double (total_number_blocks));

      ++subblock_index;
    }

    if (create_knt_mesh_parm_)
    {
      const auto knot_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      parm_block->SetBlock(subblock_index, knot_block);
      parm_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Knot mesh");
      parm_gen_->set_knot_grids(knot_block);

      this->UpdateProgress((++progress_index) / double (total_number_blocks));
    }
  } // create_parametric_mesh


  return 1;
}



void
IgatoolsParaViewReader::
set_grid_type(int arg,
              const char *const name,
              vtkGridType &type)
{

  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                << name <<  " to " << arg);

  switch (arg)
  {
    case 0:
      if (type != vtkGridType::UnstructuredQuadratic)
      {
        type = vtkGridType::UnstructuredQuadratic;
        this->Modified();
      }
      break;
    case 1:
      if (type != vtkGridType::UnstructuredLinear)
      {
        type = vtkGridType::UnstructuredLinear;
        this->Modified();
      }
      break;
    case 2:
      if (type != vtkGridType::Structured)
      {
        type = vtkGridType::Structured;
        this->Modified();
      }
      break;
    default:
      Assert(arg >= 0 && arg < 3, ExcIndexRange(arg, 0, 3));
      break;
  }
}



void
IgatoolsParaViewReader::
set_num_vis_elements(int arg1, int arg2, int arg3,
                     const char *const name,
                     const char *const mesh_type,
                     TensorSize<3> &arr)
{
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                << name << " to (" << arg1 << "," << arg2 << "," << arg3 << ")");
  if ((arr[0] != arg1) || (arr[1] != arg2) || (arr[2] != arg3))
  {
    arr[0] = arg1;
    arr[1] = arg2;
    arr[2] = arg3;
    this->Modified();
  }

  if (arg1 < 1 || arg2 < 1 || arg3 < 1)
  {
    if (arg1 < 2) arr[0] = 1;
    if (arg2 < 2) arr[1] = 1;
    if (arg3 < 2) arr[2] = 1;

    vtkWarningMacro(<< "In IgatoolsParaViewReader invalid specified number of visualization elements "
                    << "per Bezier element for the " << mesh_type << " mesh("
                    << arg1 << ", " << arg2 << ", " << arg3 << "). "
                    << "All the values must be >= 1.\n"
                    << "The number of elements was automatically set to ("
                    << arr[0] << ", " << arr[1] << ", " << arr[2] << ").\n");
  }
}



void
IgatoolsParaViewReader::
SetNumVisualizationElementsPhysicalSolid(int arg1, int arg2, int arg3)
{
  this->set_num_vis_elements(arg1, arg2, arg3,
                             "NumVisualizationElementsPhysicalSolid",
                             "physical solid",
                             n_vis_elem_phys_solid_);
}



void
IgatoolsParaViewReader::
SetNumVisualizationElementsPhysicalKnot(int arg1, int arg2, int arg3)
{
  this->set_num_vis_elements(arg1, arg2, arg3,
                             "NumVisualizationElementsPhysicalKnot",
                             "physical knot",
                             n_vis_elem_phys_knot_);
}



void
IgatoolsParaViewReader::
SetNumVisualizationElementsParametricSolid(int arg1, int arg2, int arg3)
{
  this->set_num_vis_elements(arg1, arg2, arg3,
                             "NumVisualizationElementsParametricSolid",
                             "parametric solid",
                             n_vis_elem_parm_solid_);
}



void
IgatoolsParaViewReader::
SetNumVisualizationElementsParametricKnot(int arg1, int arg2, int arg3)
{
  this->set_num_vis_elements(arg1, arg2, arg3,
                             "NumVisualizationElementsParametricKnot",
                             "parametric knot",
                             n_vis_elem_parm_knot_);
};



void
IgatoolsParaViewReader::
SetGridTypePhysicalSolid(int arg)
{
  this->set_grid_type(arg, "GridTypePhysicalSolid", phys_sol_grid_type_);
}



void
IgatoolsParaViewReader::
SetGridTypePhysicalKnot(int arg)
{
  this->set_grid_type(arg, "GridTypePhysicalKnot", phys_knt_grid_type_);
  Assert(phys_knt_grid_type_ != vtkGridType::Structured,
         ExcMessage("Knot mesh must be unstructured."));
}



void
IgatoolsParaViewReader::
SetGridTypePhysicalControl(int arg)
{
  this->set_grid_type(arg, "GridTypePhysicalControl", phys_ctr_grid_type_);
  Assert(phys_ctr_grid_type_ != vtkGridType::UnstructuredQuadratic,
         ExcMessage("Control mesh cannot be quadratic."));
}



void
IgatoolsParaViewReader::
SetGridTypeParametricSolid(int arg)
{
  this->set_grid_type(arg, "GridTypeParametricSolid", parm_sol_grid_type_);
}



void
IgatoolsParaViewReader::
SetGridTypeParametricKnot(int arg)
{
  this->set_grid_type(arg, "GridTypeParametricKnot", parm_knt_grid_type_);
  Assert(parm_knt_grid_type_ != vtkGridType::Structured,
         ExcMessage("Knot mesh must be unstructured."));
}



void
IgatoolsParaViewReader::
SetSolidMeshPhysical(bool arg)
{
  auto &name  = this->create_sol_mesh_phys_;
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                "SolidMeshPhysical"
                " to " << arg);
  if (name != arg)
  {
    name = arg;
    this->Modified();
  }
}



void
IgatoolsParaViewReader::
SetKnotMeshPhysical(bool arg)
{
  auto &name  = this->create_knt_mesh_phys_;
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                "KnotMeshPhysical"
                " to " << arg);
  if (name != arg)
  {
    name = arg;
    this->Modified();
  }
}



void
IgatoolsParaViewReader::
SetControlMeshPhysical(bool arg)
{
  auto &name  = this->create_ctr_mesh_phys_;
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                "ControlMeshPhysical"
                " to " << arg);
  if (name != arg)
  {
    name = arg;
    this->Modified();
  }
}



void
IgatoolsParaViewReader::
SetSolidMeshParametric(bool arg)
{
  auto &name  = this->create_sol_mesh_parm_;
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                "SolidMeshParametric"
                " to " << arg);
  if (name != arg)
  {
    name = arg;
    this->Modified();
  }
}



void
IgatoolsParaViewReader::
SetKnotMeshParametric(bool arg)
{
  auto &name  = this->create_knt_mesh_parm_;
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                "KnotMeshParametric"
                " to " << arg);
  if (name != arg)
  {
    name = arg;
    this->Modified();
  }
}



void
IgatoolsParaViewReader::
SetParametricMesh(bool arg)
{
  auto &name  = this->create_parametric_mesh_;
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                "create_parametric_mesh_"
                " to " << arg);
  if (name != arg)
  {
    name = arg;
    this->Modified();
  }
}



void
IgatoolsParaViewReader::
SetPhysicalMesh(bool arg)
{
  auto &name  = this->create_physical_mesh_;
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                "PhysicalMesh"
                " to " << arg);
  if (name != arg)
  {
    name = arg;
    this->Modified();
  }
}



void
IgatoolsParaViewReader::
SetFileName(const char *arg)
{
  vtkDebugMacro(<< this->GetClassName() << " (" << this <<
                "): setting FileName to  " << arg);

  if (this->file_name_ && arg && (!strcmp(this->file_name_, arg)))
    return;

  if (this->file_name_)
  {
    delete [] this->file_name_;
  }

  if (arg)
  {
    this->file_name_ = new char[strlen(arg)+1];
    strcpy(this->file_name_, arg);

    parse_file_ = true;
  }
  else
  {
    this->file_name_ = NULL;
  }

  this->Modified();
}



int
IgatoolsParaViewReader::
GetNumberOfPhysGeomArrays()
{
  Assert(phys_gen_ != nullptr, ExcNullPtr());
  return phys_gen_->get_number_grids();
}



const char *
IgatoolsParaViewReader::
GetPhysGeomArrayName(int index)
{
  Assert(phys_gen_ != nullptr, ExcNullPtr());
  const string &name = phys_gen_->get_grid_name(index);
  return name.c_str();
}



int
IgatoolsParaViewReader::
GetPhysGeomArrayStatus(const char *name)
{
  Assert(phys_gen_ != nullptr, ExcNullPtr());
  return phys_gen_->get_grid_status(string(name));
}



void
IgatoolsParaViewReader::
SetPhysGeomArrayStatus(const char *name, int en)
{
  // Note: sometimes this function is called before parsing and
  // names gotten from Previous ParaView session are parsed.
  // The if is introduced for fixing this problem.
  if (phys_gen_ != nullptr)
  {
    const auto name_str = string(name);
    if (phys_gen_->get_grid_status(name_str) != en)
    {
      phys_gen_->set_grid_status(name_str, en);
      this->Modified();
    }
  }
}



int
IgatoolsParaViewReader::
GetNumberOfParmGeomArrays()
{
  Assert(parm_gen_ != nullptr, ExcNullPtr());
  return parm_gen_->get_number_grids();
}



const char *
IgatoolsParaViewReader::
GetParmGeomArrayName(int index)
{
  Assert(parm_gen_ != nullptr, ExcNullPtr());
  const string &name = parm_gen_->get_grid_name(index);
  return name.c_str();
}



int
IgatoolsParaViewReader::
GetParmGeomArrayStatus(const char *name)
{
  Assert(parm_gen_ != nullptr, ExcNullPtr());
  return parm_gen_->get_grid_status(string(name));
}



void
IgatoolsParaViewReader::
SetParmGeomArrayStatus(const char *name, int en)
{
  // Note: sometimes this function is called before parsing and
  // names gotten from Previous ParaView session are parsed.
  // The if is introduced for fixing this problem.
  if (parm_gen_ != nullptr)
  {
    const auto name_str = string(name);
    if (parm_gen_->get_grid_status(name_str) != en)
    {
      parm_gen_->set_grid_status(name_str, en);
      this->Modified();
    }
  }
}


#include <igatools/io/reader.h>
#include <igatools/linear_algebra/epetra_vector.h>

template <int dim>
void
IgatoolsParaViewReader::
create_geometries()
{
  using std::dynamic_pointer_cast;
  using std::const_pointer_cast;

  static const int codim = dim == 1 ? 1 : 0;
  static const int space_dim = dim + codim;
  static const int range = space_dim;
  static const int rank = 1;
  static const Transformation transf = Transformation::h_grad;

//  using Fun_ = Function<dim, 0, range, rank>;
  using FunPhys_ = Function<dim, codim, range, rank>;
//  using IgFun_ = IgFunction<dim, 0, range, rank>;
  using IgFunPhys_ = IgFunction<dim, codim, range, rank>;
  using PhysSpace_ = PhysicalSpace<dim, range, rank, codim>;
//  using RefSpace_ = ReferenceSpace<dim, range, rank>;
//  using IdFun_ = IdentityFunction<dim, dim>;
  using Map = Domain<dim,codim>;
//  using Space_ = Space<dim, 0, space_dim, 1>;
  using Grid_ = Grid<dim>;

  // File names;
  const string dir_name = "/home/martinelli/tmp/paraview_plugin/";
  const string fname_0 = dir_name + "patch_0_" + std::to_string(dim) + "D.xml";
  const string fname_1 = dir_name + "patch_1_" + std::to_string(dim) + "D.xml";
  const string fname_2 = dir_name + "patch_2_" + std::to_string(dim) + "D.xml";
  const string fname_3 = dir_name + "patch_3_" + std::to_string(dim) + "D.xml";

  // Reading maps
  const shared_ptr <Map> map_0 =
    get_mapping_from_file <dim, codim> (fname_0);
  const shared_ptr <Map> map_1 =
    get_mapping_from_file <dim, codim> (fname_1);
  const shared_ptr <Map> map_2 =
    get_mapping_from_file <dim, codim> (fname_2);
  const shared_ptr <Map> map_3 =
    get_mapping_from_file <dim, codim> (fname_3);

  funcs_container_->insert_domain(map_0);
  funcs_container_->insert_domain(map_1);
  funcs_container_->insert_domain(map_2);
  funcs_container_->insert_domain(map_3);

//  Assert(false,ExcNotImplemented());

#if 0
  // Getting ig spaces.
  shared_ptr <const Space_> space_0 =
    dynamic_pointer_cast <IgFun_> (map_0)->get_ig_space();
  shared_ptr <const Space_> space_1 =
    dynamic_pointer_cast <IgFun_> (map_1)->get_ig_space();
  shared_ptr <const Space_> space_2 =
    dynamic_pointer_cast <IgFun_> (map_2)->get_ig_space();
  shared_ptr <const Space_> space_3 =
    dynamic_pointer_cast <IgFun_> (map_3)->get_ig_space();

  // Getting reference spaces.
  const shared_ptr <RefSpace_> ref_space_0 =
    const_pointer_cast <RefSpace_> (
      dynamic_pointer_cast <const RefSpace_> (space_0));
  const shared_ptr <RefSpace_> ref_space_1 =
    const_pointer_cast <RefSpace_> (
      dynamic_pointer_cast <const RefSpace_> (space_1));
  const shared_ptr <RefSpace_> ref_space_2 =
    const_pointer_cast <RefSpace_> (
      dynamic_pointer_cast <const RefSpace_> (space_2));
  const shared_ptr <RefSpace_> ref_space_3 =
    const_pointer_cast <RefSpace_> (
      dynamic_pointer_cast <const RefSpace_> (space_3));

  // Building physical spaces.
  const shared_ptr <const PhysSpace_> phys_space_0 = PhysSpace_::create(
                                                       ref_space_0, map_0);
  const shared_ptr <const PhysSpace_> phys_space_1 = PhysSpace_::create(
                                                       ref_space_1, map_1);
  const shared_ptr <const PhysSpace_> phys_space_2 = PhysSpace_::create(
                                                       ref_space_2, map_2);
  const shared_ptr <const PhysSpace_> phys_space_3 = PhysSpace_::create(
                                                       ref_space_3, map_3);

  // Getting grids;
  const shared_ptr <Grid_> grid_0 = ref_space_0->get_ptr_grid();
  const shared_ptr <Grid_> grid_1 = ref_space_1->get_ptr_grid();
  const shared_ptr <Grid_> grid_2 = ref_space_2->get_ptr_grid();
  const shared_ptr <Grid_> grid_3 = ref_space_3->get_ptr_grid();

  // Creating identity functions.
  const auto id_map_0 = IdFun_::create(grid_0);
  const auto id_map_1 = IdFun_::create(grid_1);
  const auto id_map_2 = IdFun_::create(grid_2);
  const auto id_map_3 = IdFun_::create(grid_3);

  Epetra_SerialComm comm;

  // Creating coefficients for ig physical space functions.
  auto phys_coeff_0 = EpetraTools::create_vector(*phys_space_0, "active",
                                                 comm);
  const auto dofs_dist_0 = space_0->get_ptr_const_dof_distribution();
  const Real val0 = 10.0;
  Index counter = 0;
  for (const auto &d : dofs_dist_0->get_dofs_const_view())
  {
    double val = (counter++) * 1.5 + val0;
    phys_coeff_0->ReplaceGlobalValues(1, &val, &d);
  }

  auto phys_coeff_1 = EpetraTools::create_vector(*phys_space_1, "active",
                                                 comm);
  const auto dofs_dist_1 = space_1->get_ptr_const_dof_distribution();
  const Real val1 = 11.0;
  counter = 0;
  for (const auto &d : dofs_dist_1->get_dofs_const_view())
  {
    double val = (counter++) * 1.5 + val1;
    phys_coeff_1->ReplaceGlobalValues(1, &val, &d);
  }

  auto phys_coeff_2 = EpetraTools::create_vector(*phys_space_2, "active",
                                                 comm);
  const auto dofs_dist_2 = space_2->get_ptr_const_dof_distribution();
  const Real val2 = 12.0;
  counter = 0;
  for (const auto &d : dofs_dist_2->get_dofs_const_view())
  {
    double val = (counter++) * 1.5 + val2;
    phys_coeff_2->ReplaceGlobalValues(1, &val, &d);
  }

  auto phys_coeff_3 = EpetraTools::create_vector(*phys_space_3, "active",
                                                 comm);
  const auto dofs_dist_3 = space_3->get_ptr_const_dof_distribution();
  const Real val3 = 13.0;
  counter = 0;
  for (const auto &d : dofs_dist_3->get_dofs_const_view())
  {
    double val = (counter++) * 1.5 + val3;
    phys_coeff_3->ReplaceGlobalValues(1, &val, &d);
  }

  // Creating coefficients for ig reference space functions.
  auto ref_coeff_0 = EpetraTools::create_vector(*ref_space_0, "active",
                                                comm);
  const Real ref_val0 = 20.0;
  counter = 0;
  for (const auto &d : dofs_dist_0->get_dofs_const_view())
  {
    double val = (counter++) * 1.5 + ref_val0;
    ref_coeff_0->ReplaceGlobalValues(1, &val, &d);
  }

  auto ref_coeff_1 = EpetraTools::create_vector(*ref_space_1, "active",
                                                comm);
  const Real ref_val1 = 21.0;
  counter = 0;
  for (const auto &d : dofs_dist_1->get_dofs_const_view())
  {
    double val = (counter++) * 1.5 + ref_val1;
    ref_coeff_1->ReplaceGlobalValues(1, &val, &d);
  }

  auto ref_coeff_2 = EpetraTools::create_vector(*ref_space_2, "active",
                                                comm);
  const Real ref_val2 = 22.0;
  counter = 0;
  for (const auto &d : dofs_dist_2->get_dofs_const_view())
  {
    double val = (counter++) * 1.5 + ref_val2;
    ref_coeff_2->ReplaceGlobalValues(1, &val, &d);
  }

  auto ref_coeff_3 = EpetraTools::create_vector(*ref_space_3, "active",
                                                comm);
  const Real ref_val3 = 23.0;
  counter = 0;
  for (const auto &d : dofs_dist_3->get_dofs_const_view())
  {
    double val = (counter++) * 1.5 + ref_val3;
    ref_coeff_3->ReplaceGlobalValues(1, &val, &d);
  }

  // Creating ig functions for physical spaces.
  auto ps_func_0 = dynamic_pointer_cast <FunPhys_> (
                     IgFunPhys_::create(phys_space_0, phys_coeff_0));
  auto ps_func_1 = dynamic_pointer_cast <FunPhys_> (
                     IgFunPhys_::create(phys_space_1, phys_coeff_1));
  auto ps_func_2 = dynamic_pointer_cast <FunPhys_> (
                     IgFunPhys_::create(phys_space_2, phys_coeff_2));
  auto ps_func_3 = dynamic_pointer_cast <FunPhys_> (
                     IgFunPhys_::create(phys_space_3, phys_coeff_3));

  // Creating ig functions for reference spaces.
  auto rf_func_0 = dynamic_pointer_cast <Fun_> (
                     IgFun_::create(ref_space_0, ref_coeff_0));
  auto rf_func_1 = dynamic_pointer_cast <Fun_> (
                     IgFun_::create(ref_space_1, ref_coeff_1));
  auto rf_func_2 = dynamic_pointer_cast <Fun_> (
                     IgFun_::create(ref_space_2, ref_coeff_2));
  auto rf_func_3 = dynamic_pointer_cast <Fun_> (
                     IgFun_::create(ref_space_3, ref_coeff_3));

  // Adding all the stuff to the functions container.

  // Inserting geometries.
  const string map_name_0 = "map_0_" + std::to_string(dim) + "D";
  const string map_name_1 = "map_1_" + std::to_string(dim) + "D";
  const string map_name_2 = "map_2_" + std::to_string(dim) + "D";
  const string map_name_3 = "map_3_" + std::to_string(dim) + "D";
  funcs_container_->insert_domain(
    const_pointer_cast <Map> (
      phys_space_0->get_ptr_const_map_func()),
    map_name_0);
  funcs_container_->insert_domain(
    const_pointer_cast <Map> (
      phys_space_1->get_ptr_const_map_func()),
    map_name_1);
  funcs_container_->insert_domain(
    const_pointer_cast <Map> (
      phys_space_2->get_ptr_const_map_func()),
    map_name_2);
  funcs_container_->insert_domain(
    const_pointer_cast <Map> (
      phys_space_3->get_ptr_const_map_func()),
    map_name_3);
  /*
    const string id_map_name_0 = "id_map_0_" + std::to_string(dim) + "D";
    const string id_map_name_1 = "id_map_1_" + std::to_string(dim) + "D";
    const string id_map_name_2 = "id_map_2_" + std::to_string(dim) + "D";
    const string id_map_name_3 = "id_map_3_" + std::to_string(dim) + "D";
    funcs_container_->insert_domain(id_map_0, id_map_name_0);
    funcs_container_->insert_domain(id_map_1, id_map_name_1);
    funcs_container_->insert_domain(id_map_2, id_map_name_2);
    funcs_container_->insert_domain(id_map_3, id_map_name_3);
  //*/

  // Inserting associated functions.
  const string fun_map_name_0 = "phys_func_0_" + std::to_string(dim)
                                + "D";
  const string fun_map_name_1 = "phys_func_1_" + std::to_string(dim)
                                + "D";
  const string fun_map_name_2 = "phys_func_2_" + std::to_string(dim)
                                + "D";
  const string fun_map_name_3 = "phys_func_3_" + std::to_string(dim)
                                + "D";

  funcs_container_->insert_function(
    const_pointer_cast <Map> (
      phys_space_0->get_ptr_const_map_func()),
    ps_func_0, fun_map_name_0);
  funcs_container_->insert_function(
    const_pointer_cast <Map> (
      phys_space_1->get_ptr_const_map_func()),
    ps_func_1, fun_map_name_1);
  funcs_container_->insert_function(
    const_pointer_cast <Map> (
      phys_space_2->get_ptr_const_map_func()),
    ps_func_2, fun_map_name_2);
  funcs_container_->insert_function(
    const_pointer_cast <Map> (
      phys_space_3->get_ptr_const_map_func()),
    ps_func_3, fun_map_name_3);
#endif




#if 0 // The combination dim=1, codim=1, range=2, rank=1 is not instantiated.
  const string id_fun_map_name_0 = "ref_func_0_" + std::to_string(dim) + "D";
  const string id_fun_map_name_1 = "ref_func_1_" + std::to_string(dim) + "D";
  const string id_fun_map_name_2 = "ref_func_2_" + std::to_string(dim) + "D";
  const string id_fun_map_name_3 = "ref_func_3_" + std::to_string(dim) + "D";
  funcs_container_->insert_function(id_map_0, rf_func_0, id_fun_map_name_0);
  funcs_container_->insert_function(id_map_1, rf_func_1, id_fun_map_name_1);
  funcs_container_->insert_function(id_map_2, rf_func_2, id_fun_map_name_2);
  funcs_container_->insert_function(id_map_3, rf_func_3, id_fun_map_name_3);
#endif
}
