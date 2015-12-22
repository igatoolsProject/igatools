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

#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>
#include <vtkPointSet.h>

#include <paraview_plugin/grid_generator_container.h>

#include <igatools/functions/functions_container.h>
#include <paraview_plugin/grid_information.h>

using namespace boost::fusion;

using std::string;
using std::to_string;
using std::shared_ptr;
using std::remove_reference;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

VtkIgaGridGeneratorContainer::
VtkIgaGridGeneratorContainer(const ObjContPtr_ objs_container,
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
  this->set_names();
}



auto
VtkIgaGridGeneratorContainer::
create(const ObjContPtr_ objs_container,
       const GridInfoPtr_ phys_solid_info,
       const GridInfoPtr_ phys_knot_info,
       const ControlGridInfoPtr_ phys_control_info,
       const GridInfoPtr_ parm_solid_info,
       const GridInfoPtr_ parm_knot_info) -> SelfPtr_
{
  return SelfPtr_ (new Self_(objs_container, phys_solid_info, phys_knot_info,
                             phys_control_info, parm_solid_info, parm_knot_info));
}



void
VtkIgaGridGeneratorContainer::
fill_objects_container(const ObjContPtr_ objs_container_old)
{
  using FunctionPtrs = typename ObjectsContainer::FunctionPtrs;
  using DomainPtrs = typename ObjectsContainer::DomainPtrs;
  using GridFuncPtrs = typename ObjectsContainer::GridFuncPtrs;
  using GridPtrs = typename ObjectsContainer::GridPtrs;

  // Adding all the constant and non constant functions to the new container.
  // It is checked if its domains are present in the container.
  // If not, they are included in the new container.
  FunctionPtrs valid_f_ptr_types;
  for_each(valid_f_ptr_types, [&](const auto &ptr_type)
  {
    using FunctionType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using DomainType = Domain<FunctionType::dim, FunctionType::codim>;

    // Non constant functions from the old container.
    for (const auto &id : objs_container_old->template get_object_ids<FunctionType>())
    {
        const auto func = objs_container_old->template get_object<FunctionType>(id);
        objs_container_->template insert_const_object<FunctionType> (func);

        const auto domain = func->get_domain();
        const auto dom_id = domain->get_object_id();
        if (!objs_container_old->template is_const_object_present<DomainType>(dom_id) &&
            !objs_container_old->template is_object_present<DomainType>(dom_id))
            objs_container_->template insert_const_object<DomainType> (domain);
    }

    // Constant functions from the old container.
    for (const auto &id : objs_container_old->template get_const_object_ids<FunctionType>())
    {
        const auto func = objs_container_old->template get_const_object<FunctionType>(id);

        objs_container_->template insert_const_object<FunctionType> (func);

        const auto domain = func->get_domain();
        const auto dom_id = domain->get_object_id();
        if (!objs_container_old->template is_const_object_present<DomainType>(dom_id) &&
            !objs_container_old->template is_object_present<DomainType>(dom_id))
            objs_container_->template insert_const_object<DomainType> (domain);
    }
  });


  // Adding all the constant and non constant domains to the new container.
  // It is checked if its grid functions are present in the container.
  // If not, they are included in the new container.
  DomainPtrs valid_d_ptr_types;
  for_each(valid_d_ptr_types, [&](const auto &ptr_type)
  {
    using DomainType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridFuncType = typename DomainType::GridFuncType;

    // Constant domains from the new container.
    for (const auto &id : objs_container_->template get_const_object_ids<DomainType>())
    {
        const auto domain = objs_container_old->template get_const_object<DomainType>(id);

        const auto grid_func = domain->get_grid_function();
        const auto gf_id = grid_func->get_object_id();
        if (!objs_container_old->template is_const_object_present<GridFuncType>(gf_id) &&
            !objs_container_old->template is_object_present<GridFuncType>(gf_id))
            objs_container_->template insert_const_object<GridFuncType> (grid_func);
    }

    // Non constant domains from the old container.
    for (const auto &id : objs_container_old->template get_object_ids<DomainType>())
    {
        const auto domain = objs_container_old->template get_object<DomainType>(id);
        objs_container_->template insert_const_object<DomainType> (domain);

        const auto grid_func = domain->get_grid_function();
        const auto gf_id = grid_func->get_object_id();
        if (!objs_container_old->template is_const_object_present<GridFuncType>(gf_id) &&
            !objs_container_old->template is_object_present<GridFuncType>(gf_id))
            objs_container_->template insert_const_object<GridFuncType> (grid_func);
    }

    // Constant domains from the old container.
    for (const auto &id : objs_container_old->template get_const_object_ids<DomainType>())
    {
        const auto domain = objs_container_old->template get_const_object<DomainType>(id);
        objs_container_->template insert_const_object<DomainType> (domain);

        const auto grid_func = domain->get_grid_function();
        const auto gf_id = grid_func->get_object_id();
        if (!objs_container_old->template is_const_object_present<GridFuncType>(gf_id) &&
            !objs_container_old->template is_object_present<GridFuncType>(gf_id))
            objs_container_->template insert_const_object<GridFuncType> (grid_func);
    }
  });


  // Adding all the constant and non constant grid functions to the new container.
  // It is checked if its grids are present in the container. If not,
  // they are included in the new container.
  GridFuncPtrs valid_gf_ptr_types;
  for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
  {
    using GridFuncType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridType = typename GridFuncType::GridType;

    // Constant grid functions from the new container.
    for (const auto &id : objs_container_->template get_const_object_ids<GridFuncType>())
    {
        const auto grid_func = objs_container_old->template get_const_object<GridFuncType>(id);

        const auto grid = grid_func->get_grid();
        const auto g_id = grid->get_object_id();
        if (!objs_container_old->template is_const_object_present<GridType>(g_id) &&
            !objs_container_old->template is_object_present<GridType>(g_id))
            objs_container_->template insert_const_object<GridType> (grid);
    }

    // Non constant grid functions from the old container.
    for (const auto &id : objs_container_old->template get_object_ids<GridFuncType>())
    {
        const auto grid_func = objs_container_old->template get_object<GridFuncType>(id);
        objs_container_->template insert_const_object<GridFuncType> (grid_func);

        const auto grid = grid_func->get_grid();
        const auto g_id = grid->get_object_id();
        if (!objs_container_old->template is_const_object_present<GridType>(g_id) &&
            !objs_container_old->template is_object_present<GridType>(g_id))
            objs_container_->template insert_const_object<GridType> (grid);
    }

    // Constant grid functions from the old container.
    for (const auto &id : objs_container_old->template get_const_object_ids<GridFuncType>())
    {
        const auto grid_func = objs_container_old->template get_const_object<GridFuncType>(id);
        objs_container_->template insert_const_object<GridFuncType> (grid_func);

        const auto grid = grid_func->get_grid();
        const auto g_id = grid->get_object_id();
        if (!objs_container_old->template is_const_object_present<GridType>(g_id) &&
            !objs_container_old->template is_object_present<GridType>(g_id))
            objs_container_->template insert_const_object<GridType> (grid);
    }
  });


  // Adding all the constant and non constant grids to the new container.
  GridPtrs valid_g_ptr_types;
  for_each(valid_g_ptr_types, [&](const auto &ptr_type)
  {
    using GridType = typename remove_reference<decltype(ptr_type)>::type::element_type;

    // Non-const objects.
    for (const auto &id : objs_container_old->template get_object_ids<GridType>())
        objs_container_->template insert_const_object<GridType> (objs_container_old->template get_object<GridType>(id));

    // Const objects.
    for (const auto &id : objs_container_old->template get_const_object_ids<GridType>())
        objs_container_->template insert_const_object<GridType> (objs_container_old->template get_const_object<GridType>(id));
  });
}



void
VtkIgaGridGeneratorContainer::
set_names()
{
  using FunctionPtrs = typename ObjectsContainer::FunctionPtrs;
  using DomainPtrs = typename ObjectsContainer::DomainPtrs;
  using GridFuncPtrs = typename ObjectsContainer::GridFuncPtrs;
  using GridPtrs = typename ObjectsContainer::GridPtrs;

  FunctionPtrs valid_f_ptr_types;
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

  DomainPtrs valid_d_ptr_types;
  for_each(valid_d_ptr_types, [&](const auto &ptr_type)
  {
    using DomainType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridFuncType = typename DomainType::GridFuncType;
    using GridType = typename GridFuncType::GridType;

    for (const auto &id : objs_container_->template get_const_object_ids<DomainType>())
    {
        const auto domain = const_pointer_cast<DomainType>(
                objs_container_->template get_const_object<DomainType>(id));

        if (domain->get_name() == "")
            domain->set_name("Domain Id=" + to_string(domain->get_object_id()));

        const auto grid_func = const_pointer_cast<GridFuncType>(domain->get_grid_function());
        const auto grid = const_pointer_cast<GridType>(grid_func->get_grid());

        if (grid_func->get_name() == "")
            grid_func->set_name("GridFunction of the Domain \"" + domain->get_name() + "\"");

        if (grid->get_name() == "")
            grid->set_name("Grid of the Domain \"" + domain->get_name() + "\"");
    }
  });

  GridFuncPtrs valid_gf_ptr_types;
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

  GridPtrs valid_g_ptr_types;
  for_each(valid_g_ptr_types, [&](const auto &ptr_type)
  {
    using GridType = typename remove_reference<decltype(ptr_type)>::type::element_type;

    for (const auto &id : objs_container_->template get_const_object_ids<GridType>())
    {
        const auto grid = const_pointer_cast<GridType>(
                objs_container_->template get_const_object<GridType>(id));

        if (grid->get_name() == "")
            grid->set_name("Grid Id=" + to_string(grid->get_object_id()));
    }
  });
}



void
VtkIgaGridGeneratorContainer::
update(const GridInfoPtr_ phys_solid_info,
       const GridInfoPtr_ phys_knot_info,
       const ControlGridInfoPtr_ phys_control_info,
       const GridInfoPtr_ parm_solid_info,
       const GridInfoPtr_ parm_knot_info)
{
  AssertThrow (false, ExcNotImplemented());

  const bool phys_solid_updated = phys_solid_info_->update(phys_solid_info);
  const bool phys_knot_updated  = phys_knot_info_->update(phys_knot_info);
  const bool phys_control_updated  = phys_control_info_->update(phys_control_info);

  // Updating physical generators.
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gen_pair : gen_vec)
          gen_pair.second->update(phys_solid_updated, phys_knot_updated,
                                  phys_control_updated);
  });

  const bool parm_solid_updated = phys_solid_info_->update(parm_solid_info);
  const bool parm_knot_updated  = phys_knot_info_->update(parm_knot_info);

  // Updating parametric generators.
  boost::fusion::for_each(parm_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gen_pair : gen_vec)
          gen_pair.second->update(parm_solid_updated, parm_knot_updated,
                                  false);
  });

}



Size
VtkIgaGridGeneratorContainer::
get_number_physical_grids() const
{
  Size counter = 0;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_vec)
  {
      counter += gen_vec.size();
  });
  return counter;
}



Size
VtkIgaGridGeneratorContainer::
get_number_parametric_grids() const
{
  Size counter = 0;
  boost::fusion::for_each(parm_generators_, [&](const auto &gen_vec)
  {
      counter += gen_vec.size();
  });
  return counter;
}



Size
VtkIgaGridGeneratorContainer::
get_number_active_physical_grids() const
{
  Size counter = 0;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gp : gen_vec)
      {
          if (gp.first.is_active())
              ++counter;
      }
  });
  return counter;
}



Size
VtkIgaGridGeneratorContainer::
get_number_active_parametric_grids() const
{
  Size counter = 0;
  boost::fusion::for_each(parm_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gp : gen_vec)
      {
          if (gp.first.is_active())
              ++counter;
      }
  });
  return counter;
}



Size
VtkIgaGridGeneratorContainer::
get_number_active_ig_grids() const
{
  Size counter = 0;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gp : gen_vec)
      {
          if (gp.first.is_active() && gp.first.is_ig_grid_func())
              ++counter;
      }
  });
  return counter;
}



string
VtkIgaGridGeneratorContainer::
get_physical_grid_name(const Index &id) const
{
  string name = "";
  bool found = false;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_vec)
  {
      if (found)
          return;

      for (const auto &gp : gen_vec)
      {
          if (gp.first.get_id() == id)
          {
              found = true;
              name = gp.first.get_name();
              return;
          }
      }
  });

  Assert (found, ExcMessage("Not present id."));

  return name;
}



string
VtkIgaGridGeneratorContainer::
get_parametric_grid_name(const Index &id) const
{
  string name = "";
  bool found = false;
  boost::fusion::for_each(parm_generators_, [&](const auto &gen_vec)
  {
      if (found)
          return;

      for (const auto &gp : gen_vec)
      {
          if (gp.first.get_id() == id)
          {
              found = true;
              name = gp.first.get_name();
              return;
          }
      }
  });

  Assert (found, ExcMessage("Not present id."));

  return name;
}



bool
VtkIgaGridGeneratorContainer::
get_physical_grid_status(const std::string &name) const
{
  bool status = false;
  bool found = false;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_vec)
  {
      if (found)
          return;

      for (const auto &gp : gen_vec)
      {
          if (gp.first.get_name() == name)
          {
              found = true;
              status = gp.first.is_active();
              return;
          }
      }
  });

  Assert (found, ExcMessage("Not present name."));

  return status;
}



bool
VtkIgaGridGeneratorContainer::
get_parametric_grid_status(const std::string &name) const
{
  bool status = false;
  bool found = false;
  boost::fusion::for_each(parm_generators_, [&](const auto &gen_vec)
  {
      if (found)
          return;

      for (const auto &gp : gen_vec)
      {
          if (gp.first.get_name() == name)
          {
              found = true;
              status = gp.first.is_active();
              return;
          }
      }
  });

  Assert (found, ExcMessage("Not present name."));

  return status;
}



void
VtkIgaGridGeneratorContainer::
set_physical_grid_status(const std::string &name, const bool status)
{
  bool found = false;
  boost::fusion::for_each(phys_generators_, [&](auto &gen_vec)
  {
      if (found)
          return;

      for (auto &gp : gen_vec)
      {
          GridInfo &grid_info = gp.first;
          if (grid_info.get_name() == name)
          {
              found = true;
              grid_info.set_status(status);
          }
      }
  });

  Assert (found, ExcMessage("Not present name."));
}



void
VtkIgaGridGeneratorContainer::
set_parametric_grid_status(const std::string &name, const bool status)
{
  bool found = false;
  boost::fusion::for_each(parm_generators_, [&](auto &gen_vec)
  {
      if (found)
          return;

      for (auto &gp : gen_vec)
      {
          GridInfo &grid_info = gp.first;
          if (grid_info.get_name() == name)
          {
              found = true;
              grid_info.set_status(status);
          }
      }
  });

  Assert (found, ExcMessage("Not present name."));
}



void
VtkIgaGridGeneratorContainer::
set_physical_solid_grids(vtkMultiBlockDataSet *const mb)
{
  const auto active_grids = this->get_number_active_physical_grids();
  Assert(active_grids > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_grids);

  unsigned int block_index = 0;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gp : gen_vec)
      {
          const auto &grid_info = gp.first;

          if (grid_info.is_active())
          {
              const auto &name = grid_info.get_name();

              mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
              mb->SetBlock(block_index, gp.second->get_solid_grid());
              ++block_index;
          }
      }
  });
}



void
VtkIgaGridGeneratorContainer::
set_physical_knot_grids(vtkMultiBlockDataSet *const mb)
{
  const auto active_grids = this->get_number_active_physical_grids();
  Assert(active_grids > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_grids);

  unsigned int block_index = 0;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gp : gen_vec)
      {
          const auto &grid_info = gp.first;

          if (grid_info.is_active())
          {
              const auto &name = grid_info.get_name();

              mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
              mb->SetBlock(block_index, gp.second->get_knot_grid());
              ++block_index;
          }
      }
  });
}



void
VtkIgaGridGeneratorContainer::
set_physical_control_grids(vtkMultiBlockDataSet *const mb)
{
  const auto active_grids = this->get_number_active_ig_grids();
  Assert(active_grids > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_grids);

  unsigned int block_index = 0;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gp : gen_vec)
      {
          const auto &grid_info = gp.first;

          if (grid_info.is_active() && grid_info.is_ig_grid_func())
          {
              const auto &name = grid_info.get_name();

              mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
              mb->SetBlock(block_index, gp.second->get_control_grid());
              ++block_index;
          }
      }
  });
}



void
VtkIgaGridGeneratorContainer::
set_parametric_solid_grids(vtkMultiBlockDataSet *const mb)
{
  const auto active_grids = this->get_number_active_parametric_grids();
  Assert(active_grids > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_grids);

  unsigned int block_index = 0;
  boost::fusion::for_each(parm_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gp : gen_vec)
      {
          const auto &grid_info = gp.first;

          if (grid_info.is_active())
          {
              const auto &name = grid_info.get_name();

              mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
              mb->SetBlock(block_index, gp.second->get_solid_grid());
              ++block_index;
          }
      }
  });
}



void
VtkIgaGridGeneratorContainer::
set_parametric_knot_grids(vtkMultiBlockDataSet *const mb)
{
  const auto active_grids = this->get_number_active_parametric_grids();
  Assert(active_grids > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_grids);

  unsigned int block_index = 0;
  boost::fusion::for_each(parm_generators_, [&](const auto &gen_vec)
  {
      for (const auto &gp : gen_vec)
      {
          const auto &grid_info = gp.first;

          if (grid_info.is_active())
          {
              const auto &name = grid_info.get_name();

              mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
              mb->SetBlock(block_index, gp.second->get_knot_grid());
              ++block_index;
          }
      }
  });
}



#if 0

void
VtkIgaGridGeneratorBu::
fill_generators()
{
  boost::fusion::for_each(generators_, [&](const auto &generators_vec)
  {
      for (const auto &generator : generators_vec)
      {
//          for (const auto &domain : domains)
//          {
//              this->insert_generator(domain.get_ptr_const_data());
//              // Is it a physical mapping?
//              /*
//        if (std::dynamic_pointer_cast<IdFun>(map_fun) == nullptr)
//          this->insert_generator<dim, codim> (map_fun, name);
//        //*/
//          } // endl loop on map_funs with a given pair <dim,codim>
      }
  });
}



template<int dim, int codim>
void
VtkIgaGridGeneratorBu::
insert_generator(const DomainPtr_<dim, codim> domain)
{
#if 0
  if (this->is_physical())
  {
    Assert(domain != nullptr, ExcNullPtr());


    using IgGridFun = IgGridFunction<dim,dim + codim>;

    const bool is_ig_mapping = std::dynamic_pointer_cast<const IgGridFun>(domain->get_grid_function()) != nullptr;
    AssertThrow(is_ig_mapping, ExcMessage("Not an Isogeometric mapping!"));

    // Inserting the new generators indices in the table.
    generators_numbering_.push_back(
      make_tuple(domain->get_object_id(), domain->get_name(), true, is_ig_mapping));

    auto &data_same_dim_codim = this->template get_data_dim_codim<dim, codim>();
    auto &generators = data_same_dim_codim.get_generators();

    const Index map_id = generators_numbering_.size() - 1;

    Assert(generators.find(map_id) == generators.end(),
           ExcMessage("Key already introduced."));

    generators[map_id] = VtkIgaGridGeneratorPhys<dim, codim>::create
                         (domain, solid_info_, knot_info_, control_info_, objs_container_);

  }
  else
  {
    // Inserting the new generators indices in the table.
    generators_numbering_.push_back(
      make_tuple(domain->get_object_id(), domain->get_name(), true, false));

    auto &data_same_dim_codim = this->template get_data_dim_codim<dim, codim>();
    auto &generators = data_same_dim_codim.get_generators();

    const Index map_id = generators_numbering_.size() - 1;

    Assert(generators.find(map_id) == generators.end(),
           ExcMessage("Key already introduced."));

    generators[map_id] = VtkIgaGridGeneratorParm<dim, codim>::create
                         (domain, solid_info_, knot_info_, objs_container_);
  }
#endif
}



void
VtkIgaGridGeneratorContParm::
fill_generators()
{
  boost::fusion::for_each(generators_, [&](const auto &generators_vec)
  {
      for (const auto &generator : generators_vec)
      {
//          for (const auto &domain : domains)
//          {
//              this->insert_generator(domain.get_ptr_const_data());
//              /*
//        // Is it an identity mapping?
//        if (std::dynamic_pointer_cast<IdFun>(map_fun) != nullptr)
//          this->insert_generator<dim, codim> (map_fun, name);
//        //*/
//          } // endl loop on map_funs with a given pair <dim,codim>
      }
  });
}


template<int dim, int codim>
void
VtkIgaGridGeneratorContParm::
insert_generator(const DomainPtr_<dim, codim> domain)
{
#if 0
  Assert(domain != nullptr, ExcNullPtr());

  // Inserting the new generators indices in the table.
  generators_numbering_.push_back(
    make_tuple(domain->get_object_id(), domain->get_name(), true, false));

  auto &data_same_dim_codim = this->template get_data_dim_codim<dim, codim>();
  auto &generators = data_same_dim_codim.get_generators();

  const Index map_id = generators_numbering_.size() - 1;

  Assert(generators.find(map_id) == generators.end(),
         ExcMessage("Key already introduced."));

  generators[map_id] = VtkIgaGridGeneratorParm<dim, codim>::create
                       (domain, solid_info_, knot_info_, objs_container_);
#endif
}



#endif

IGA_NAMESPACE_CLOSE
