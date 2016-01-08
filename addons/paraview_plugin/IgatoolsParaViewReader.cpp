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

#include <igatools/io/xml_document.h>
#include <igatools/io/objects_container_xml_reader.h>
#include <igatools/base/objects_container.h>
#include <paraview_plugin/vtk_iga_grid_container.h>
#include <paraview_plugin/vtk_iga_grid_information.h>

#include <sys/stat.h>


using std::get;
using std::string;
using std::shared_ptr;
using namespace iga;

vtkStandardNewMacro(IgatoolsParaViewReader);


#ifndef SERIALIZATION
#ifndef XML_IO
static_assert(true, "Neither serialization nor XML capabilities are active.");
#endif
#endif

IgatoolsParaViewReader::IgatoolsParaViewReader()
  :
  n_vis_elem_phys_solid_(1),
  n_vis_elem_parm_solid_(1),
  n_vis_elem_phys_knot_(1),
  n_vis_elem_parm_knot_(1),
  grid_gen_(GridGenPtr_())
{
#ifndef NDEBUG
  this->DebugOn();
#else
  this->DebugOff();
#endif

  this->SetNumberOfInputPorts(0); // No vtk input, this is not a filter.
  this->SetNumberOfOutputPorts(1); // Just one output, a multi block.
}



const char *
IgatoolsParaViewReader::
GetClassNameInternal() const
{
    return "IgatoolsParaViewReader";
}



int
IgatoolsParaViewReader::
IsTypeOf(const char *type)
{
    if ( !strcmp("IgatoolsParaViewReader", type) )
    {
        return 1;
    }
    return vtkMultiBlockDataSetAlgorithm::IsTypeOf(type);
}



int
IgatoolsParaViewReader::
IsA(const char *type)
{
    return this->IgatoolsParaViewReader::IsTypeOf(type);
}



IgatoolsParaViewReader *
IgatoolsParaViewReader::
SafeDownCast(vtkObjectBase *o)
{
    if ( o && o->IsA("IgatoolsParaViewReader") )
    {
        return static_cast<IgatoolsParaViewReader *>(o); \
    }
    return NULL;
}



IgatoolsParaViewReader *
IgatoolsParaViewReader::
NewInstance() const
{
    return IgatoolsParaViewReader::SafeDownCast(this->NewInstanceInternal());
}



vtkObjectBase *
IgatoolsParaViewReader::
NewInstanceInternal() const
{
    return IgatoolsParaViewReader::New();
}



int IgatoolsParaViewReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *output_vec)
{
  vtkInformation *info = output_vec->GetInformationObject(0);

  if (!vtkMultiBlockDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT())))
    return 0;

  // If the file is not parse, it is parsed now.
  if (parse_file_)
    return this->parse_file();

  return 1;
}



int IgatoolsParaViewReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *output_vec)
{
  vtkInformation *info = output_vec->GetInformationObject(0);
  vtkMultiBlockDataSet *output =
    vtkMultiBlockDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));

  this->SetProgressText("Generating igatools geometries.");

  this->UpdateProgress(0.0);

  this->update_grid_info();
  return this->create_grids(output);
}



void IgatoolsParaViewReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: "
     << (this->file_name_ ? this->file_name_ : "(none)") << "\n";
}



int
IgatoolsParaViewReader::
CanReadFile(const char *name)
{
  // TODO this can be done if XML_IO is active
  try
  {
    XMLDocument::check_file(name);
    return 1;
  }
  catch (ExceptionBase &exc)
  {
    std::ostringstream stream;
    exc.print_exc_data(stream);
    vtkErrorMacro(<< stream.str());
    return 0;
  }
  catch (...)
  {
    vtkErrorMacro(<< "Impossible to read IGA file " << name << ".");
    return 0;
  }
}



int
IgatoolsParaViewReader::
parse_file()
{
  if (this->CanReadFile(file_name_) == 0)
    return 0;

  const auto file_name_str = string(file_name_);

  grid_gen_.reset();
  objs_container_.reset();

  // Physical solid grid.
  const auto phys_sol = VtkGridInformation::create
                        (n_vis_elem_phys_solid_, phys_sol_grid_type_);

  // Physical knot grid.
  const auto phys_knt = VtkGridInformation::create
                        (n_vis_elem_phys_knot_, phys_knt_grid_type_);

  // Physical control grid.
  const auto phys_ctr = VtkControlGridInformation::create
                        (phys_ctr_grid_type_ == VtkGridType::Structured);

  // Parametric solid grid.
  const auto parm_sol = VtkGridInformation::create
                        (n_vis_elem_parm_solid_, parm_sol_grid_type_);

  // Parametric knot grid.
  const auto parm_knt = VtkGridInformation::create
                        (n_vis_elem_parm_knot_, parm_knt_grid_type_);

  try
  {
    this->SetProgressText("Parsing igatools file.");

#ifdef XML_IO

    // TODO: before parsing the hole file (that can be big),
    // it is checked if the file has the expected structure (at least
    // the header) for knowing if it is an XML human readable or
    // a serialized file.

    // Check here if the file is of type XML human readable
    const bool xml_human_readable = true;
    if (xml_human_readable)
    {
      objs_container_ = ObjectsContainerXMLReader::parse_const(file_name_str);
      AssertThrow(!objs_container_->is_void(),
                  ExcMessage("No objects defined in the input file: "
                             + file_name_str + "."));
      parse_file_ = false;
    }

#endif

#ifdef SERIALIZATION
    if (parse_file_)
    {
      ObjectsContainer container_new;
      {
        std::ifstream xml_istream(file_name_str);
        IArchive xml_in(xml_istream);
        xml_in >> container_new;
      }
      objs_container_ = std::make_shared<ObjectsContainer>(container_new);

      AssertThrow(!objs_container_->is_void(),
                  ExcMessage("No objects defined in the input file or "
                             "serialization file not properly defined."
                             " File name: " + file_name_str + "."));
      parse_file_ = false;
    }
#endif

    grid_gen_ = VtkIgaGridContainer::create
                (objs_container_, phys_sol, phys_knt, phys_ctr,
                 parm_sol, parm_knt);

    return 1;
  }
  catch (std::exception &e)
  {
    vtkErrorMacro(<< e.what());

    objs_container_ = ObjectsContainer::create();
    grid_gen_ = VtkIgaGridContainer::create
                (objs_container_, phys_sol, phys_knt, phys_ctr,
                 parm_sol, parm_knt);

    return 0;
  }
  catch (...)
  {
    vtkErrorMacro(<< "An exception occurred when parsing file "
                  << file_name_str << ".");

    objs_container_ = ObjectsContainer::create();
    grid_gen_ = VtkIgaGridContainer::create
                (objs_container_, phys_sol, phys_knt, phys_ctr,
                 parm_sol, parm_knt);

    return 0;
  }
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
                        (phys_ctr_grid_type_ == VtkGridType::Structured);

  // Parametric solid grid.
  const auto parm_sol = VtkGridInformation::create
                        (n_vis_elem_parm_solid_, parm_sol_grid_type_);


  // Parametric knot grid.
  const auto parm_knt = VtkGridInformation::create
                        (n_vis_elem_parm_knot_, parm_knt_grid_type_);

  grid_gen_->update(phys_sol, phys_knt, phys_ctr, parm_sol, parm_knt);
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

  const auto num_active_phys = grid_gen_->get_number_active_physical_grids();
  const auto num_active_parm = grid_gen_->get_number_active_parametric_grids();

  bool new_create_physical_mesh = create_physical_mesh_;
  if (create_physical_mesh_ && (num_phys_blocks == 0 || num_active_phys == 0))
  {
    if (num_phys_blocks == 0)
    {
      vtkWarningMacro(<< "Physical geometries set active, but no grid type "
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
      vtkWarningMacro(<< "Parametric geometries set active, but no grid type "
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
  {
    vtkWarningMacro(<< "Neither physical nor parametric geometries are "
                    "active. No output produced.");

    return 1;
  }


  // Creating blocks for the physical and parametric geometries.
  mb->SetNumberOfBlocks(num_blocks);
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
      grid_gen_->set_physical_solid_grids(solid_block);

      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));

      ++subblock_index;
    }

    if (create_knt_mesh_phys_)
    {
      const auto knot_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      phys_block->SetBlock(subblock_index, knot_block);
      phys_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Knot mesh");
      grid_gen_->set_physical_knot_grids(knot_block);

      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));

      ++subblock_index;
    }

    if (create_ctr_mesh_phys_)
    {
      const auto control_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      phys_block->SetBlock(subblock_index, control_block);
      phys_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Control mesh");
      grid_gen_->set_physical_control_grids(control_block);

      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));
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
      grid_gen_->set_parametric_solid_grids(solid_block);

      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));

      ++subblock_index;
    }

    if (create_knt_mesh_parm_)
    {
      const auto knot_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      parm_block->SetBlock(subblock_index, knot_block);
      parm_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Knot mesh");
      grid_gen_->set_parametric_knot_grids(knot_block);

      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));
    }
  } // create_parametric_mesh


  return 1;
}



void
IgatoolsParaViewReader::
set_grid_type(int arg,
              const char *const name,
              VtkGridType &type)
{

  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting "
                << name <<  " to " << arg);

  switch (arg)
  {
    case 0:
      if (type != VtkGridType::UnstructuredQuadratic)
      {
        type = VtkGridType::UnstructuredQuadratic;
        this->Modified();
      }
      break;
    case 1:
      if (type != VtkGridType::UnstructuredLinear)
      {
        type = VtkGridType::UnstructuredLinear;
        this->Modified();
      }
      break;
    case 2:
      if (type != VtkGridType::Structured)
      {
        type = VtkGridType::Structured;
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

    vtkWarningMacro(<< "In IgatoolsParaViewReader invalid specified "
                    << "number of visualization elements per Bezier "
                    << "for the " << mesh_type << " mesh(" << arg1
                    << ", " << arg2 << ", " << arg3 << "). All the values"
                    << " must be >= 1.\nThe number of elements was "
                    << "automatically set to (" << arr[0] << ", "
                    << arr[1] << ", " << arr[2] << ").\n");
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
  Assert(phys_knt_grid_type_ != VtkGridType::Structured,
         ExcMessage("Knot mesh must be unstructured."));
}



void
IgatoolsParaViewReader::
SetGridTypePhysicalControl(int arg)
{
  this->set_grid_type(arg, "GridTypePhysicalControl", phys_ctr_grid_type_);
  Assert(phys_ctr_grid_type_ != VtkGridType::UnstructuredQuadratic,
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
  Assert(parm_knt_grid_type_ != VtkGridType::Structured,
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
  if (grid_gen_ == nullptr)
    return 0;
  return grid_gen_->get_number_physical_grids();
}



const char *
IgatoolsParaViewReader::
GetPhysGeomArrayName(int index)
{
  Assert(grid_gen_ != nullptr, ExcNullPtr());
  const char *name = grid_gen_->get_physical_grid_name(index);
  return name;
}



int
IgatoolsParaViewReader::
GetPhysGeomArrayStatus(const char *name)
{
  Assert(grid_gen_ != nullptr, ExcNullPtr());
  return grid_gen_->get_physical_grid_status(string(name));
}



void
IgatoolsParaViewReader::
SetPhysGeomArrayStatus(const char *name, int enable)
{
  // Note: sometimes this function is called before parsing and
  // names gotten from Previous ParaView session are parsed.
  // The if is introduced for fixing this problem.
  if (grid_gen_ != nullptr)
  {
    const auto name_str = string(name);
    if (grid_gen_->get_physical_grid_status(name_str) != enable)
    {
      grid_gen_->set_physical_grid_status(name_str, enable);
      this->Modified();
    }
  }
}



int
IgatoolsParaViewReader::
GetNumberOfParmGeomArrays()
{
  if (grid_gen_ == nullptr)
    return 0;
  return grid_gen_->get_number_parametric_grids();
}



const char *
IgatoolsParaViewReader::
GetParmGeomArrayName(int index)
{
  Assert(grid_gen_ != nullptr, ExcNullPtr());
  const char *name = grid_gen_->get_parametric_grid_name(index);
  return name;
}



int
IgatoolsParaViewReader::
GetParmGeomArrayStatus(const char *name)
{
  Assert(grid_gen_ != nullptr, ExcNullPtr());
  return grid_gen_->get_parametric_grid_status(string(name));
}



void
IgatoolsParaViewReader::
SetParmGeomArrayStatus(const char *name, int enable)
{
  // Note: sometimes this function is called before parsing and
  // names gotten from Previous ParaView session are parsed.
  // The if is introduced for fixing this problem.
  if (grid_gen_ != nullptr)
  {
    const auto name_str = string(name);
    if (grid_gen_->get_parametric_grid_status(name_str) != enable)
    {
      grid_gen_->set_parametric_grid_status(name_str, enable);
      this->Modified();
    }
  }
}
