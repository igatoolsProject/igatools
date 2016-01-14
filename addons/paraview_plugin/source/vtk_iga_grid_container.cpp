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

#include <paraview_plugin/vtk_iga_grid.h>

#include <paraview_plugin/vtk_iga_grid_container.h>
#include <paraview_plugin/vtk_iga_grid_information.h>

#include <igatools/functions/grid_function_lib.h>
#include <igatools/io/objects_container_xml_reader.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>
#include <vtkPointSet.h>

#include <sys/stat.h>


using namespace boost::fusion;

using std::string;
using std::to_string;
using std::shared_ptr;
using std::remove_reference;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

namespace paraview_plugin
{

VtkIgaGridContainer::
VtkIgaGridContainer(const ObjContPtr_ objs_container,
                    const GridInfoPtr_ phys_solid_info,
                    const GridInfoPtr_ phys_knot_info,
                    const ControlGridInfoPtr_ phys_control_info,
                    const GridInfoPtr_ parm_solid_info,
                    const GridInfoPtr_ parm_knot_info)
  :
  objs_container_(ObjectsContainer::create()),
  phys_solid_info_(phys_solid_info),
  phys_knot_info_(phys_knot_info),
  phys_control_info_(phys_control_info),
  parm_solid_info_(parm_solid_info),
  parm_knot_info_(parm_knot_info)
{
  Assert(objs_container != nullptr, ExcNullPtr());
  Assert(objs_container_ != nullptr, ExcNullPtr());
  Assert(phys_solid_info_ != nullptr, ExcNullPtr());
  Assert(phys_knot_info_ != nullptr, ExcNullPtr());
  Assert(phys_control_info_ != nullptr, ExcNullPtr());
  Assert(parm_solid_info_ != nullptr, ExcNullPtr());
  Assert(parm_knot_info_ != nullptr, ExcNullPtr());

  this->fill_objects_container(objs_container);
  this->set_container_names();
  this->build_vtk_iga_grids();

  // Checking if there is any valid VTK IGA grid.
  auto lambda_func_not_empty = [](const auto &gen_pair)
  {
    return !gen_pair.second.empty();
  };

  if (!boost::fusion::any(phys_grids_, lambda_func_not_empty) &&
      !boost::fusion::any(parm_grids_, lambda_func_not_empty))
      // Both lists of VTK IGA grids are empty.
      throw ExcVtkError("There are not valid domains or grids.");
}



VtkIgaGridContainer::
VtkIgaGridContainer()
  :
  objs_container_(ObjectsContainer::create()),
  phys_solid_info_(VtkGridInformation::create_void()),
  phys_knot_info_(VtkGridInformation::create_void()),
  phys_control_info_(VtkControlGridInformation::create_void()),
  parm_solid_info_(VtkGridInformation::create_void()),
  parm_knot_info_(VtkGridInformation::create_void())
{
  Assert(objs_container_ != nullptr, ExcNullPtr());
  Assert(phys_solid_info_ != nullptr, ExcNullPtr());
  Assert(phys_knot_info_ != nullptr, ExcNullPtr());
  Assert(phys_control_info_ != nullptr, ExcNullPtr());
  Assert(parm_solid_info_ != nullptr, ExcNullPtr());
  Assert(parm_knot_info_ != nullptr, ExcNullPtr());

  this->build_vtk_iga_grids();
}


auto
VtkIgaGridContainer::
create(const std::string &file_name,
       const NumCells_   &n_cells_phs_sol,
       const VtkGridType &grid_type_phs_sol,
       const NumCells_   &n_cells_phs_knt,
       const VtkGridType &grid_type_phs_knt,
       const VtkGridType &grid_type_phs_ctr,
       const NumCells_   &n_cells_prm_sol,
       const VtkGridType &grid_type_prm_sol,
       const NumCells_   &n_cells_prm_knt,
       const VtkGridType &grid_type_prm_knt) -> SelfPtr_
{
  // Physical solid grid.
  const auto phys_sol = VtkGridInformation::create
                        (n_cells_phs_sol, grid_type_phs_sol);

  // Physical knot grid.
  const auto phys_knt = VtkGridInformation::create
                        (n_cells_phs_knt, grid_type_phs_knt);

  // Physical control grid.
  const auto phys_ctr = VtkControlGridInformation::create
                        (grid_type_phs_ctr == VtkGridType::Structured);

  // Parametric solid grid.
  const auto parm_sol = VtkGridInformation::create
                        (n_cells_prm_sol, grid_type_prm_sol);

  // Parametric knot grid.
  const auto parm_knt = VtkGridInformation::create
                        (n_cells_prm_knt, grid_type_prm_knt);

  const auto objs_container = parse_objects_container(file_name);

  return SelfPtr_(new Self_(objs_container, phys_sol, phys_knt, phys_ctr,
                            parm_sol, parm_knt));

}



auto
VtkIgaGridContainer::
parse_objects_container(const string &file_name) ->
ObjContPtr_
{
    Self_::check_file(file_name);

    ObjContPtr_ objs_container;

    if (Self_::is_file_binary(file_name))
    {
#ifdef SERIALIZATION
        ObjectsContainer container_new;
        {
            std::ifstream xml_istream(file_name);
            IArchive xml_in(xml_istream);
            xml_in >> container_new;
        }
        objs_container = std::make_shared<ObjectsContainer>(container_new);
#endif
    }
    else
    {
#ifdef XML_IO
        objs_container = ObjectsContainerXMLReader::parse_const(file_name);
#endif
    }

    Assert (objs_container != nullptr, ExcNullPtr());
    return objs_container;
}



void
VtkIgaGridContainer::
check_file(const std::string &file_name)
{
  // Checking if the file exists.
  errno = 0;
  struct stat buffer;
  if (stat(file_name.c_str(), &buffer) == -1) // == 0 ok; == -1 error
  {
    if (errno == ENOENT)   // errno declared by include file errno.h
      throw ExcVtkError("Path file does not exist, or path is an empty string.");
    else if (errno == ENOTDIR)
      throw ExcVtkError("A component of the path is not a directory.");
    else if (errno == ELOOP)
      throw ExcVtkError("Too many symbolic links encountered while traversing the path.");
    else if (errno == EACCES)
      throw ExcVtkError("Permission denied.");
    else if (errno == ENAMETOOLONG)
      throw ExcVtkError("File can not be read");
    else
      throw ExcVtkError("An unknown problem was encountered.");
  }

    if (Self_::is_file_binary(file_name))
    {
#ifndef SERIALIZATION
        throw ExcVtkError("Impossible to parse binary format file."
                " Igatools serialization of binary files was not activated."
                " Currently this plugin is only capable of parsing XML files.");
#endif
    }
    else
    {
#ifndef XML_IO
        throw ExcVtkError("Impossible to parse an ascii format file."
                " Igatools XML ascii parser was not activated."
                " Currently this plugin is only capable of parsing "
                "igatools serialized files in binary format.");
#endif
    }

}



bool
VtkIgaGridContainer::
is_file_binary(const std::string &file_name)
{
    int c;
    std::ifstream file(file_name);
    while((c = file.get()) != EOF && c <= 127);
    file.close();
    if(c == EOF)
        return false;
    else
        return true;
}



SafeSTLVector<string>
VtkIgaGridContainer::
get_invalid_dimension_objects(const ObjContPtr_ objs_container)
{
  SafeSTLVector<string> invalid_objects;

  InvalidGridPtrs_ invalid_g_ptr_types;
  for_each(invalid_g_ptr_types, [&](const auto &ptr_type)
           {
              using ObjectType = typename remove_reference<decltype(ptr_type)>::type::element_type;
              for (const auto &id : objs_container->template get_const_object_ids<ObjectType>())
              {
                  const auto obj = objs_container->template get_const_object<ObjectType>(id);
                  invalid_objects.push_back(
                          "Grid<" + to_string(ObjectType::dim) + ">, "
                          "Name: " + obj->get_name() + ", "
                          "ObjectId: " + to_string(id));
              }
           });


  InvalidDomainPtrs_ invalid_d_ptr_types;
  for_each(invalid_d_ptr_types, [&](const auto &ptr_type)
           {
              using ObjectType = typename remove_reference<decltype(ptr_type)>::type::element_type;
              for (const auto &id : objs_container->template get_const_object_ids<ObjectType>())
              {
                  const auto obj = objs_container->template get_const_object<ObjectType>(id);
                  invalid_objects.push_back(
                          "Domain<" + to_string(ObjectType::dim) + ", " +
                          to_string(ObjectType::space_dim - ObjectType::dim) + ">, "
                          "Name: " + obj->get_name() + ", "
                          "ObjectId: " + to_string(id));
              }
           });

  InvalidGridFuncPtrs_ invalid_gf_ptr_types;
  for_each(invalid_gf_ptr_types, [&](const auto &ptr_type)
           {
              using ObjectType = typename remove_reference<decltype(ptr_type)>::type::element_type;
              for (const auto &id : objs_container->template get_const_object_ids<ObjectType>())
              {
                  const auto obj = objs_container->template get_const_object<ObjectType>(id);
                  invalid_objects.push_back(
                          "GridFunction<" + to_string(ObjectType::dim) + ", " +
                          to_string(ObjectType::range) + ">, "
                          "Name: " + obj->get_name() + ", "
                          "ObjectId: " + to_string(id));
              }
           });

  InvalidFunctionPtrs_ invalid_f_ptr_types;
  for_each(invalid_f_ptr_types, [&](const auto &ptr_type)
           {
              using ObjectType = typename remove_reference<decltype(ptr_type)>::type::element_type;
              for (const auto &id : objs_container->template get_const_object_ids<ObjectType>())
              {
                  const auto obj = objs_container->template get_const_object<ObjectType>(id);
                  invalid_objects.push_back(
                          "Function<" + to_string(ObjectType::dim) + ", "
                                      + to_string(ObjectType::codim) + ", "
                                      + to_string(ObjectType::range) + ", "
                                      + to_string(ObjectType::rank) + ">, "
                          "Name: " + obj->get_name() + ", "
                          "ObjectId: " + to_string(id));
              }
           });

  return invalid_objects;
}



auto
VtkIgaGridContainer::
create_void() -> SelfPtr_
{
  return SelfPtr_(new Self_());
}



void
VtkIgaGridContainer::
create_multiblock_grid(const bool phys_mesh,
                       const bool sol_phys_mesh,
                       const bool knt_phys_mesh,
                       const bool ctr_phys_mesh,
                       const bool prm_mesh,
                       const bool sol_prm_mesh,
                       const bool knt_prm_mesh,
                       vtkMultiBlockDataSet *const mb) const
{
  unsigned int num_blocks = phys_mesh + prm_mesh;

  if (num_blocks == 0)
  {
      VtkIgaWarningMacro("Neither physical nor parametric geometries are "
                         "active. No output produced.");
  }

  const unsigned int num_phys_blocks = sol_phys_mesh + knt_phys_mesh + ctr_phys_mesh;

  const unsigned int num_parm_blocks = sol_prm_mesh + knt_prm_mesh;

  // TO improve
  const auto num_active_phys = Self_::get_number_active_domains(phys_grids_);
  const auto num_active_parm = Self_::get_number_active_domains(parm_grids_);

  bool new_create_physical_mesh = phys_mesh;
  if (phys_mesh && (num_phys_blocks == 0 || num_active_phys == 0))
  {
      if (num_phys_blocks == 0)
      {
          VtkIgaWarningMacro("Physical geometries set active, but no grid type "
                  "(solid, knot, control) has been selected");
      }
      else
      {
          VtkIgaWarningMacro("Physical geometries set active, but no "
                  "geometry set active from the list.");
      }

    --num_blocks;

    new_create_physical_mesh = false;
  }


  bool new_create_parametric_mesh = prm_mesh;
  if (prm_mesh && (num_parm_blocks == 0 || num_active_parm == 0))
  {
      if (num_parm_blocks == 0)
      {
          VtkIgaWarningMacro("Parametric geometries set active, but no grid type "
                  "(solid, knot) has been selected");
      }
      else
      {
          VtkIgaWarningMacro("Parametric geometries set active, but no "
                  "geometry set active from the list.");
      }

    --num_blocks;

    new_create_parametric_mesh = false;
  }

  if (num_blocks == 0)
  {
      VtkIgaWarningMacro("Neither physical nor parametric geometries are "
              "active. No output produced.");
  }


  // Creating blocks for the physical and parametric geometries.
  mb->SetNumberOfBlocks(num_blocks);
  for (unsigned int i = 0; i < num_blocks; ++i)
    mb->SetBlock(i, vtkSmartPointer<vtkMultiBlockDataSet>::New());

  Size total_number_blocks = 0;
  if (new_create_physical_mesh)
  {
    if (sol_phys_mesh)
      ++total_number_blocks;
    if (knt_phys_mesh)
      ++total_number_blocks;
    if (ctr_phys_mesh)
      ++total_number_blocks;
  }
  if (new_create_parametric_mesh)
  {
    if (sol_prm_mesh)
      ++total_number_blocks;
    if (knt_prm_mesh)
      ++total_number_blocks;
  }


//  Index progress_index = 0;

  unsigned int block_index = 0;
  if (new_create_physical_mesh)
  {
    mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), "Physical domains");

    vtkMultiBlockDataSet *const phys_block =
      vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index));

    phys_block->SetNumberOfBlocks(num_phys_blocks);

    Index subblock_index = 0;

    if (sol_phys_mesh)
    {
      const auto solid_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      phys_block->SetBlock(subblock_index, solid_block);
      phys_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Solid meshes");
      Self_::set_solid_domains(phys_grids_, solid_block);

//      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));

      ++subblock_index;
    }

    if (knt_phys_mesh)
    {
      const auto knot_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      phys_block->SetBlock(subblock_index, knot_block);
      phys_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Knot meshes");
      Self_::set_knot_domains(phys_grids_, knot_block);

//      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));

      ++subblock_index;
    }

    if (ctr_phys_mesh)
    {
      const auto control_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      phys_block->SetBlock(subblock_index, control_block);
      phys_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Control polygon meshes");
      Self_::set_control_domains(phys_grids_, control_block);

//      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));
    }

    ++block_index;
  } // create_physical_mesh


  if (new_create_parametric_mesh)
  {
    mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), "Parametric domains");

    vtkMultiBlockDataSet *const parm_block =
      vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index));

    parm_block->SetNumberOfBlocks(num_parm_blocks);

    Index subblock_index = 0;

    if (sol_prm_mesh)
    {
      const auto solid_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      parm_block->SetBlock(subblock_index, solid_block);
      parm_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Solid meshes");
      Self_::set_solid_domains(parm_grids_, solid_block);

//      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));

      ++subblock_index;
    }

    if (knt_prm_mesh)
    {
      const auto knot_block = vtkSmartPointer <vtkMultiBlockDataSet>::New();
      parm_block->SetBlock(subblock_index, knot_block);
      parm_block->GetMetaData(subblock_index)->Set(vtkCompositeDataSet::NAME(),
                                                   "Knot meshes");
      Self_::set_knot_domains(parm_grids_, knot_block);

//      this->UpdateProgress(double (++progress_index) / double (total_number_blocks));
    }
  } // create_parametric_mesh

}



void
VtkIgaGridContainer::
update(const NumCells_   &n_cells_phs_sol,
       const VtkGridType &grid_type_phs_sol,
       const NumCells_   &n_cells_phs_knt,
       const VtkGridType &grid_type_phs_knt,
       const VtkGridType &grid_type_phs_ctr,
       const NumCells_   &n_cells_prm_sol,
       const VtkGridType &grid_type_prm_sol,
       const NumCells_   &n_cells_prm_knt,
       const VtkGridType &grid_type_prm_knt)
{
  // Physical solid grid.
  const auto phys_solid_info = VtkGridInformation::create
                        (n_cells_phs_sol, grid_type_phs_sol);

  // Physical knot grid.
  const auto phys_knot_info = VtkGridInformation::create
                        (n_cells_phs_knt, grid_type_phs_knt);

  // Physical control grid.
  const auto phys_control_info = VtkControlGridInformation::create
                        (grid_type_phs_ctr == VtkGridType::Structured);

  // Parametric solid grid.
  const auto parm_solid_info = VtkGridInformation::create
                        (n_cells_prm_sol, grid_type_prm_sol);

  // Parametric knot grid.
  const auto parm_knot_info = VtkGridInformation::create
                        (n_cells_prm_knt, grid_type_prm_knt);

  const bool phys_solid_updated    = phys_solid_info_->update(phys_solid_info);
  const bool phys_knot_updated     = phys_knot_info_->update(phys_knot_info);
  const bool phys_control_updated  = phys_control_info_->update(phys_control_info);
  const bool parm_solid_updated    = parm_solid_info_->update(parm_solid_info);
  const bool parm_knot_updated     = parm_knot_info_->update(parm_knot_info);
  const bool parm_control_updated  = false;

  // Updating physical VTK IGA grids.
  boost::fusion::for_each(phys_grids_, [&](const auto &gen_pair)
  {
    for (const auto gen : gen_pair.second)
      gen->update(phys_solid_updated, phys_knot_updated,
                  phys_control_updated);
  });


  // Updating parametric VTK IGA grids.
  boost::fusion::for_each(parm_grids_, [&](const auto &gen_pair)
  {
    for (const auto &gen : gen_pair.second)
      gen->update(parm_solid_updated, parm_knot_updated,
                  parm_control_updated);
  });
}



void
VtkIgaGridContainer::
fill_objects_container(const ObjContPtr_ objs_container_old)
{
    // Adding missing object into the container.
    objs_container_old->fill_not_inserted_dependencies();

    // Checking if there is any object with invalid dimension in the
    // container.
    const auto invalid_dim_obj = Self_::get_invalid_dimension_objects(objs_container_old);
    if (invalid_dim_obj.size() > 0)
    {
        string warning_message = "The following objects have dimensions that the "
                "plugin cannot currently visualize:\n";
        for (const auto &obj : invalid_dim_obj)
            warning_message += obj + "\n";

        VtkIgaWarningMacro(<< warning_message);
    }

    // Adding all the present objects in the container to a new container,
    // all of them as constant.
    auto lambda_insert_as_const = [&](const auto &ptr_type)
    {
        using Type = typename remove_reference<decltype(ptr_type)>::type::element_type;
        for (const auto &id : objs_container_old->template get_object_ids<Type>())
            objs_container_->template insert_const_object<Type>(
                    objs_container_old->template get_object<Type>(id));

        for (const auto &id : objs_container_old->template get_const_object_ids<Type>())
            objs_container_->template insert_const_object<Type>(
                    objs_container_old->template get_const_object<Type>(id));
    };

    GridPtrs_ valid_g_ptr_types;
    for_each(valid_g_ptr_types, lambda_insert_as_const);

    GridFuncPtrs_ valid_gf_ptr_types;
    for_each(valid_gf_ptr_types, lambda_insert_as_const);

    DomainPtrs_ valid_d_ptr_types;
    for_each(valid_d_ptr_types, lambda_insert_as_const);

    FunctionPtrs_ valid_f_ptr_types;
    for_each(valid_f_ptr_types, lambda_insert_as_const);
}



void
VtkIgaGridContainer::
set_container_names()
{
  // Dummy objects, just for iterate dynamically over all the types.
  FunctionPtrs_ valid_f_ptr_types;
  DomainPtrs_ valid_d_ptr_types;
  GridFuncPtrs_ valid_gf_ptr_types;
  GridPtrs_ valid_g_ptr_types;


  for_each(valid_f_ptr_types, [&](const auto &ptr_type)
  {
    using FunctionType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using DomainType = Domain<FunctionType::dim, FunctionType::codim>;
    using GridFuncType = typename DomainType::GridFuncType;
    using GridType = typename GridFuncType::GridType;

    for (const auto &id : objs_container_->template get_const_object_ids<FunctionType>())
    {
      const auto func = const_pointer_cast<FunctionType>(
                          objs_container_->template get_const_object<FunctionType>(id));

      if (func->get_name() == "")
        func->set_name("Function Id=" + to_string(func->get_object_id()));

      const auto domain = const_pointer_cast<DomainType>(func->get_domain());
      const auto grid_func = const_pointer_cast<GridFuncType>(domain->get_grid_function());
      const auto grid = const_pointer_cast<GridType>(grid_func->get_grid());

      if (domain->get_name() == "")
      {
        if (grid_func->get_name() == "")
          domain->set_name("Domain of Function \"" + func->get_name() + "\"");
        else
          domain->set_name("Domain with GridFunction \"" + grid_func->get_name() + "\"");
      }

      if (grid_func->get_name() == "")
        grid_func->set_name("GridFunction of the Domain \"" + domain->get_name() + "\"");

      if (grid->get_name() == "")
        grid->set_name("Grid of the Domain \"" + domain->get_name() + "\"");
    }
  });


  std::set<string> domain_names;

  for_each(valid_d_ptr_types, [&](const auto &ptr_type)
  {
    using DomainType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridFuncType = typename DomainType::GridFuncType;
    using GridType = typename GridFuncType::GridType;
    using ValidFuncPtrs_ = ValidFuncsForDomain<DomainType>;

    ValidFuncPtrs_ valid_fd_ptr_types;

    for (const auto &id : objs_container_->template get_const_object_ids<DomainType>())
    {
      const auto domain = const_pointer_cast<DomainType>(objs_container_->template get_const_object<DomainType>(id));

      if (domain->get_name() == "")
        domain->set_name("Domain Id=" + to_string(domain->get_object_id()));

      if (domain_names.find(domain->get_name()) != domain_names.cend())
      {
        Index i = 1;
        string name = domain->get_name() + " (" + to_string(++i) + ")";
        while (domain_names.find(name) != domain_names.cend())
        {
          name = domain->get_name() + " (" + to_string(++i) + ")";
        }
        domain->set_name(name);
      }

      domain_names.insert(domain->get_name());

      const auto grid_func = const_pointer_cast<GridFuncType>(domain->get_grid_function());
      const auto grid = const_pointer_cast<GridType>(grid_func->get_grid());

      if (grid_func->get_name() == "")
        grid_func->set_name("GridFunction of the Domain \"" + domain->get_name() + "\"");

      if (grid->get_name() == "")
        grid->set_name("Grid of the Domain \"" + domain->get_name() + "\"");

      // Checking for repeated function names associated to this domain.
      std::set<string> function_names;

      for_each(valid_fd_ptr_types, [&](const auto &f_ptr_type)
      {
        using FunctionType = typename remove_reference<decltype(f_ptr_type)>::type::element_type;

        for (const auto &f_id : objs_container_->template get_const_object_ids<FunctionType>())
        {
          const auto function = const_pointer_cast<FunctionType>(objs_container_->template get_const_object<FunctionType>(f_id));
          const auto &name = function->get_name();

          if (function->get_domain() == domain)
          {
            if (function_names.find(name) != function_names.cend())
            {
              Index i = 1;
              string new_name = name + " (" + to_string(++i) + ")";
              while (function_names.find(new_name) != function_names.cend())
              {
                new_name = name + " (" + to_string(++i) + ")";
              }

              function->set_name(new_name);
            }

            function_names.insert(function->get_name());
          }
        }

      });

    }
  });

  for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
  {
    using GridFuncType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridType = typename GridFuncType::GridType;

    for (const auto &id : objs_container_->template get_const_object_ids<GridFuncType>())
    {
      const auto grid_func = const_pointer_cast<GridFuncType>(
                               objs_container_->template get_const_object<GridFuncType>(id));

      if (grid_func->get_name() == "")
        grid_func->set_name("GridFunction Id=" + to_string(grid_func->get_object_id()));

      const auto grid = const_pointer_cast<GridType>(grid_func->get_grid());
      if (grid->get_name() == "")
        grid->set_name("Grid of the GridFunction \"" + grid_func->get_name() + "\"");
    }
  });


  std::set<string> grid_names;

  for_each(valid_g_ptr_types, [&](const auto &ptr_type)
  {
    using GridType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using ValidGridFuncPtrs_ = ValidGridFuncsForDomain<Domain<GridType::dim, 0>>;

    ValidGridFuncPtrs_ valid_gfd_ptr_types;

    for (const auto &id : objs_container_->template get_const_object_ids<GridType>())
    {
      const auto grid = const_pointer_cast<GridType>(
                          objs_container_->template get_const_object<GridType>(id));

      if (grid->get_name() == "")
        grid->set_name("Grid Id=" + to_string(grid->get_object_id()));

      if (grid_names.find(grid->get_name()) != grid_names.cend())
      {
        Index i = 1;
        string name = grid->get_name() + " (" + to_string(++i) + ")";
        while (grid_names.find(name) != grid_names.cend())
        {
          name = grid->get_name() + " (" + to_string(++i) + ")";
        }
        grid->set_name(name);
      }

      grid_names.insert(grid->get_name());

      // Checking for repeated function names associated to this grid.
      std::set<string> grid_func_names;

      for_each(valid_gfd_ptr_types, [&](const auto &f_ptr_type)
      {
        using GridFuncType = typename remove_reference<decltype(f_ptr_type)>::type::element_type;

        for (const auto &gf_id : objs_container_->template get_const_object_ids<GridFuncType>())
        {
          const auto grid_func = const_pointer_cast<GridFuncType>(objs_container_->template get_const_object<GridFuncType>(gf_id));
          const auto &name = grid_func->get_name();

          if (grid_func->get_grid() == grid)
          {
            if (grid_func_names.find(name) != grid_func_names.cend())
            {
              Index i = 1;
              string new_name = name + " (" + to_string(++i) + ")";
              while (grid_func_names.find(new_name) != grid_func_names.cend())
              {
                new_name = name + " (" + to_string(++i) + ")";
              }

              grid_func->set_name(new_name);
            }

            grid_func_names.insert(grid_func->get_name());
          }
        }

      });
    }
  });
}



void
VtkIgaGridContainer::
build_vtk_iga_grids()
{
  // Building physical VTK IGA grids.
  Index grid_id = 0;
  bool is_physical = true;

  DomainPtrs_ valid_d_ptr_types;
  for_each(valid_d_ptr_types, [&](const auto &ptr_type)
  {
    using DomainType = typename remove_reference<decltype(ptr_type)>::type::element_type;

    for (const auto &id : objs_container_->template get_const_object_ids<DomainType>())
    {
      const auto domain = objs_container_->template get_const_object<DomainType>(id);

      const bool is_active = true;

      const auto gg = VtkIgaGrid<DomainType>::
                      create(domain,  grid_id, phys_solid_info_,
                             phys_knot_info_, phys_control_info_,
                             objs_container_, is_active, is_physical);

      at_key<DomainType>(phys_grids_).push_back(gg);

      ++grid_id;
    }
  });


  // Building parametric VTK IGA grids.
  grid_id = 0;
  is_physical = false;

  GridPtrs_ valid_g_ptr_types;
  for_each(valid_g_ptr_types, [&](const auto &ptr_type)
  {
    using GridType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    static const int dim = GridType::dim;
    using DomainType = Domain<dim, 0>;

    for (const auto &id : objs_container_->template get_const_object_ids<GridType>())
    {
      const auto grid = objs_container_->template get_const_object<GridType>(id);

      const bool is_active = true;

      // Creating an identity domain for the grid.
      const auto id_func = grid_functions::IdentityGridFunction<dim>::const_create(grid);
      const auto domain = DomainType::const_create(id_func);
      const_pointer_cast<DomainType>(domain)->set_name(grid->get_name());

      const auto gg = VtkIgaGrid<DomainType>::
                      create(domain, grid_id, parm_solid_info_,
                             parm_knot_info_, phys_control_info_,
                             objs_container_, is_active, is_physical);

      at_key<DomainType>(parm_grids_).push_back(gg);

      ++grid_id;
    }
  });
}



Size
VtkIgaGridContainer::
get_number_physical_domains() const
{
  return Self_::get_number_domains(phys_grids_);
}



Size
VtkIgaGridContainer::
get_number_parametric_domains() const
{
  return Self_::get_number_domains(parm_grids_);
}



const char *
VtkIgaGridContainer::
get_physical_domain_name(const Index &id) const
{
  return Self_::get_domain_name(phys_grids_, id);
}



const char *
VtkIgaGridContainer::
get_parametric_domain_name(const Index &id) const
{
  return Self_::get_domain_name(parm_grids_, id);
}



bool
VtkIgaGridContainer::
get_physical_domain_status(const std::string &name) const
{
  return Self_::get_domain_status(phys_grids_, name);
}



bool
VtkIgaGridContainer::
get_parametric_domain_status(const std::string &name) const
{
  return Self_::get_domain_status(parm_grids_, name);
}



void
VtkIgaGridContainer::
set_physical_domain_status(const std::string &name, const bool status)
{
  Self_::set_domain_status(phys_grids_, name, status);
}



void
VtkIgaGridContainer::
set_parametric_domain_status(const std::string &name, const bool status)
{
  Self_::set_domain_status(parm_grids_, name, status);
}



Size
VtkIgaGridContainer::
get_number_domains(const VtkIgaGridsContainer_ grid_container)
{
  Size counter = 0;
  boost::fusion::for_each(grid_container, [&](const auto &gen_pair)
  {
    counter += gen_pair.second.size();
  });
  return counter;
}



Size
VtkIgaGridContainer::
get_number_active_domains(const VtkIgaGridsContainer_ grid_container)
{
  Size counter = 0;
  boost::fusion::for_each(grid_container, [&](const auto &gen_pair)
  {
    for (const auto gen : gen_pair.second)
      counter += gen->is_active();
  });
  return counter;
}



const char *
VtkIgaGridContainer::
get_domain_name(const VtkIgaGridsContainer_ grid_container,
              const Index &id)
{
  const char *name = "";
  bool found = false;
  boost::fusion::for_each(grid_container, [&](const auto &gen_pair)
  {
    if (found)
      return;

    for (const auto gen : gen_pair.second)
    {
      if (gen->get_id() == id)
      {
        found = true;
        const string &name_str = gen->get_name();
        name = name_str.c_str();
        return;
      }
    }
  });

  Assert(found, ExcMessage("Not present id."));

  return name;

}



bool
VtkIgaGridContainer::
get_domain_status(const VtkIgaGridsContainer_ grid_container,
                const std::string &name)
{
  bool status = false;
  bool found = false;
  boost::fusion::for_each(grid_container, [&](const auto &gen_pair)
  {
    if (found)
      return;

    for (const auto gen : gen_pair.second)
    {
      if (gen->get_name() == name)
      {
        found = true;
        status = gen->is_active();
        return;
      }
    }
  });

  Assert(found, ExcMessage("Not present name."));

  return status;

}



void
VtkIgaGridContainer::
set_domain_status(const VtkIgaGridsContainer_ grid_container,
                const std::string &name,
                const bool status)
{
  bool found = false;
  boost::fusion::for_each(grid_container, [&](auto &gen_pair)
  {
    if (found)
      return;

    for (auto gen : gen_pair.second)
    {
      if (gen->get_name() == name)
      {
        found = true;
        gen->set_status(status);
      }
    }
  });

  Assert(found, ExcMessage("Not present name."));
}



void
VtkIgaGridContainer::
set_solid_domains(const VtkIgaGridsContainer_ grid_container,
                  vtkMultiBlockDataSet *const mb)
{
  const auto active_domains = get_number_active_domains(grid_container);
  Assert(active_domains > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_domains);

  unsigned int block_index = 0;
  boost::fusion::for_each(grid_container, [&](const auto &gen_pair)
  {
    for (const auto gen : gen_pair.second)
    {
      if (gen->is_active())
      {
        const auto &name = gen->get_name();

        mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
        mb->SetBlock(block_index, gen->get_solid_grid());
        ++block_index;
      }
    }
  });

}



void
VtkIgaGridContainer::
set_knot_domains(const VtkIgaGridsContainer_ grid_container,
               vtkMultiBlockDataSet *const mb)
{
  const auto active_domains = get_number_active_domains(grid_container);
  Assert(active_domains > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_domains);

  unsigned int block_index = 0;
  boost::fusion::for_each(grid_container, [&](const auto &gen_pair)
  {
    for (const auto gen : gen_pair.second)
    {
      if (gen->is_active())
      {
        const auto &name = gen->get_name();

        mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
        mb->SetBlock(block_index, gen->get_knot_grid());
        ++block_index;
      }
    }
  });
}



void
VtkIgaGridContainer::
set_control_domains(const VtkIgaGridsContainer_ grid_container,
                  vtkMultiBlockDataSet *const mb)
{
  // Getting number of active ig grids.
  Size active_domains = 0;
  boost::fusion::for_each(grid_container, [&](const auto &gen_pair)
  {
    for (const auto gen : gen_pair.second)
      active_domains += (gen->is_active() && gen->is_ig_grid_func());
  });

  Assert(active_domains > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_domains);

  unsigned int block_index = 0;
  boost::fusion::for_each(grid_container, [&](const auto &gen_pair)
  {
    for (const auto gen : gen_pair.second)
    {
      if (gen->is_active() && gen->is_ig_grid_func())
      {
        const auto &name = gen->get_name();

        mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
        mb->SetBlock(block_index, gen->get_control_grid());
        ++block_index;
      }
    }
  });
}

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE
