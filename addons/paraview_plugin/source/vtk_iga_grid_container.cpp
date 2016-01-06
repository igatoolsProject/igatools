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
#include <igatools/geometry/grid_function_lib.h>

#include <paraview_plugin/vtk_iga_grid_container.h>
#include <paraview_plugin/vtk_iga_grid_information.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>
#include <vtkPointSet.h>


using namespace boost::fusion;

using std::string;
using std::to_string;
using std::shared_ptr;
using std::remove_reference;
using std::const_pointer_cast;

IGA_NAMESPACE_OPEN

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
  this->set_names();
  this->build_generators();
}



auto
VtkIgaGridContainer::
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
VtkIgaGridContainer::
fill_objects_container(const ObjContPtr_ objs_container_old)
{
  // Adding all the constant and non constant functions to the new container.
  // It is checked if its domains are present in the container.
  // If not, they are included in the new container.

  // TODO: to document further here.

  // All the present functions are inserted into the new container.
  // The domain of every function is obtained. If it is not included
  // in the old container, it is inserted directly into the new container.

  FunctionPtrs_ valid_f_ptr_types;
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


  // Adding all the domains to the new container.
  // It is checked if its grid functions are present into the container.
  // If not, they are inserted into the new container.
  DomainPtrs_ valid_d_ptr_types;
  for_each(valid_d_ptr_types, [&](const auto &ptr_type)
  {
    using DomainType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridFuncType = typename DomainType::GridFuncType;

    // Inserting grid functions from the previously inserted domains.
    for (const auto &id : objs_container_->template get_const_object_ids<DomainType>())
    {
        const auto domain = objs_container_->template get_const_object<DomainType>(id);

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


  // Adding all the grid functions to the new container.
  // It is checked if its grids are present in the container. If not,
  // they are included in the new container.
  GridFuncPtrs_ valid_gf_ptr_types;
  for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
  {
    using GridFuncType = typename remove_reference<decltype(ptr_type)>::type::element_type;
    using GridType = typename GridFuncType::GridType;

    // Inserting grids from the previously inserted grid functions.
    for (const auto &id : objs_container_->template get_const_object_ids<GridFuncType>())
    {
        const auto grid_func = objs_container_->template get_const_object<GridFuncType>(id);

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


  // Adding all the grids to the new container.
  GridPtrs_ valid_g_ptr_types;
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
VtkIgaGridContainer::
set_names()
{
  FunctionPtrs_ valid_f_ptr_types;
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

  DomainPtrs_ valid_d_ptr_types;
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

  GridFuncPtrs_ valid_gf_ptr_types;
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

  GridPtrs_ valid_g_ptr_types;
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
VtkIgaGridContainer::
build_generators()
{
  // Building physical generators.
  Index domain_id = 0;
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
                create(domain,  domain_id, phys_solid_info_,
                       phys_knot_info_, phys_control_info_,
                       objs_container_, is_active, is_physical);

        at_key<DomainType>(phys_generators_).push_back(gg);

        ++domain_id;
    }
  });


  // Building parametric generators.
  domain_id = 0;
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
                create(domain, domain_id, parm_solid_info_,
                       parm_knot_info_, phys_control_info_,
                       objs_container_, is_active, is_physical);

        at_key<DomainType>(parm_generators_).push_back(gg);

        ++domain_id;
    }
  });
}



void
VtkIgaGridContainer::
update(const GridInfoPtr_ phys_solid_info,
       const GridInfoPtr_ phys_knot_info,
       const ControlGridInfoPtr_ phys_control_info,
       const GridInfoPtr_ parm_solid_info,
       const GridInfoPtr_ parm_knot_info)
{
  const bool phys_solid_updated    = phys_solid_info_->update(phys_solid_info);
  const bool phys_knot_updated     = phys_knot_info_->update(phys_knot_info);
  const bool phys_control_updated  = phys_control_info_->update(phys_control_info);
  const bool parm_solid_updated    = parm_solid_info_->update(parm_solid_info);
  const bool parm_knot_updated     = parm_knot_info_->update(parm_knot_info);
  const bool parm_control_updated  = false;

  // Updating physical generators.
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_pair)
  {
      for (const auto gen : gen_pair.second)
          gen->update(phys_solid_updated, phys_knot_updated,
                      phys_control_updated);
  });


  // Updating parametric generators.
  boost::fusion::for_each(parm_generators_, [&](const auto &gen_pair)
  {
      for (const auto &gen : gen_pair.second)
          gen->update(parm_solid_updated, parm_knot_updated,
                      parm_control_updated);
  });
}



Size
VtkIgaGridContainer::
get_number_physical_grids() const
{
  return Self_::get_number_grids(phys_generators_);
}



Size
VtkIgaGridContainer::
get_number_parametric_grids() const
{
  return Self_::get_number_grids(parm_generators_);
}



Size
VtkIgaGridContainer::
get_number_active_physical_grids() const
{
  return Self_::get_number_active_grids(phys_generators_);
}



Size
VtkIgaGridContainer::
get_number_active_parametric_grids() const
{
  return Self_::get_number_active_grids(parm_generators_);
}



Size
VtkIgaGridContainer::
get_number_active_ig_grids() const
{
  Size counter = 0;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_pair)
  {
      for (const auto gen : gen_pair.second)
          counter += (gen->is_active() && gen->is_ig_grid_func());
  });
  return counter;
}



const char *
VtkIgaGridContainer::
get_physical_grid_name(const Index &id) const
{
    return Self_::get_grid_name(phys_generators_, id);
}



const char *
VtkIgaGridContainer::
get_parametric_grid_name(const Index &id) const
{
    return Self_::get_grid_name(parm_generators_, id);
}



bool
VtkIgaGridContainer::
get_physical_grid_status(const std::string &name) const
{
    return Self_::get_grid_status(phys_generators_, name);
}



bool
VtkIgaGridContainer::
get_parametric_grid_status(const std::string &name) const
{
    return Self_::get_grid_status(parm_generators_, name);
}



void
VtkIgaGridContainer::
set_physical_grid_status(const std::string &name, const bool status)
{
    Self_::set_grid_status(phys_generators_, name, status);
}



void
VtkIgaGridContainer::
set_parametric_grid_status(const std::string &name, const bool status)
{
    Self_::set_grid_status(parm_generators_, name, status);
}



void
VtkIgaGridContainer::
set_physical_solid_grids(vtkMultiBlockDataSet *const mb)
{
  Self_::set_solid_grids(phys_generators_, mb);
}



void
VtkIgaGridContainer::
set_physical_knot_grids(vtkMultiBlockDataSet *const mb)
{
  Self_::set_knot_grids(phys_generators_, mb);
}



void
VtkIgaGridContainer::
set_physical_control_grids(vtkMultiBlockDataSet *const mb)
{
  const auto active_grids = this->get_number_active_ig_grids();
  Assert(active_grids > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_grids);

  unsigned int block_index = 0;
  boost::fusion::for_each(phys_generators_, [&](const auto &gen_pair)
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



void
VtkIgaGridContainer::
set_parametric_solid_grids(vtkMultiBlockDataSet *const mb)
{
  Self_::set_solid_grids(parm_generators_, mb);
}



void
VtkIgaGridContainer::
set_parametric_knot_grids(vtkMultiBlockDataSet *const mb)
{
  Self_::set_knot_grids(parm_generators_, mb);
}



Size
VtkIgaGridContainer::
get_number_grids(const GridGensContainer_ generators)
{
  Size counter = 0;
  boost::fusion::for_each(generators, [&](const auto &gen_pair)
  {
      counter += gen_pair.second.size();
  });
  return counter;
}



Size
VtkIgaGridContainer::
get_number_active_grids(const GridGensContainer_ generators)
{
  Size counter = 0;
  boost::fusion::for_each(generators, [&](const auto &gen_pair)
  {
      for (const auto gen : gen_pair.second)
          counter += gen->is_active();
  });
  return counter;

}



const char *
VtkIgaGridContainer::
get_grid_name(const GridGensContainer_ generators,
              const Index &id)
{
  const char *name = "";
  bool found = false;
  boost::fusion::for_each(generators, [&](const auto &gen_pair)
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

  Assert (found, ExcMessage("Not present id."));

  return name;

}



bool
VtkIgaGridContainer::
get_grid_status(const GridGensContainer_ generators,
                const std::string &name)
{
  bool status = false;
  bool found = false;
  boost::fusion::for_each(generators, [&](const auto &gen_pair)
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

  Assert (found, ExcMessage("Not present name."));

  return status;

}



void
VtkIgaGridContainer::
set_grid_status(const GridGensContainer_ generators,
                const std::string &name,
                const bool status)
{
  bool found = false;
  boost::fusion::for_each(generators, [&](auto &gen_pair)
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

  Assert (found, ExcMessage("Not present name."));
}



void
VtkIgaGridContainer::
set_solid_grids(const GridGensContainer_ generators,
                vtkMultiBlockDataSet *const mb)
{
  const auto active_grids = get_number_active_grids(generators);
  Assert(active_grids > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_grids);

  unsigned int block_index = 0;
  boost::fusion::for_each(generators, [&](const auto &gen_pair)
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
set_knot_grids(const GridGensContainer_ generators,
               vtkMultiBlockDataSet *const mb)
{
  const auto active_grids = get_number_active_grids(generators);
  Assert(active_grids > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(active_grids);

  unsigned int block_index = 0;
  boost::fusion::for_each(generators, [&](const auto &gen_pair)
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

IGA_NAMESPACE_CLOSE
