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
using std::shared_ptr;

IGA_NAMESPACE_OPEN

VtkIgaGridGeneratorContainer::
VtkIgaGridGeneratorContainer(const ObjContPtr_ objs_container,
                             const GridInfoPtr_ phys_solid_info,
                             const GridInfoPtr_ phys_knot_info,
                             const ControlGridInfoPtr_ phys_control_info,
                             const GridInfoPtr_ parm_solid_info,
                             const GridInfoPtr_ parm_knot_info)
  :
  objs_container_(objs_container),
  phys_solid_info_(phys_solid_info),
  phys_knot_info_(phys_knot_info),
  phys_control_info_(phys_control_info),
  parm_solid_info_(parm_solid_info),
  parm_knot_info_(parm_knot_info)
{
  Assert(objs_container_ != nullptr, ExcNullPtr());
  Assert(phys_solid_info_ != nullptr, ExcNullPtr());
  Assert(phys_knot_info_ != nullptr, ExcNullPtr());
  Assert(phys_control_info_ != nullptr, ExcNullPtr());
  Assert(parm_solid_info_ != nullptr, ExcNullPtr());
  Assert(parm_knot_info_ != nullptr, ExcNullPtr());

  this->fill_generators();
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
fill_generators()
{
    AssertThrow (false, ExcNotImplemented());
  // Fill container properly
  //
  // 1 - For every function, insert the domain
  // 2 - For every domain, insert its grid.
  // 2 - For every grid function, insert its grid.
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
