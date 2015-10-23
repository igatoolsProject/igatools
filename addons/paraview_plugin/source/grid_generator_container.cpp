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




VtkIgaGridGeneratorContBase::
VtkIgaGridGeneratorContBase(const FunContPtr_ funcs_container,
                            const GridInfoPtr_ solid_info,
                            const GridInfoPtr_ knot_info,
                            const ControlGridInfoPtr_ control_info)
  :
  funcs_container_(funcs_container),
  solid_info_(solid_info),
  knot_info_(knot_info),
  control_info_(control_info)
{
  Assert(funcs_container_ != nullptr, ExcNullPtr());
  Assert(solid_info_ != nullptr, ExcNullPtr());
  Assert(knot_info_ != nullptr, ExcNullPtr());

}






template<int dim, int codim>
auto
VtkIgaGridGeneratorContBase::
get_data_dim_codim() const ->
const VtkGridGenContSameCodim_<dim, codim> &
{
  return this->template get_data_dim<dim>().template get_data_codim<codim>();
}



template<int dim, int codim>
auto
VtkIgaGridGeneratorContBase::
get_data_dim_codim() ->
VtkGridGenContSameCodim_<dim, codim> &
{
  return this->template get_data_dim<dim>().template get_data_codim<codim>();
}



template <int dim>
template<int codim>
auto
VtkIgaGridGeneratorContBase::
VtkGridGeneratorContBaseSameDim<dim>::
VtkGridGeneratorContBaseSameDimCodim<codim>::
get_generators() const -> const GenMap_<dim, codim> &
{
  return grid_generators_;
}



template <int dim>
template<int codim>
auto
VtkIgaGridGeneratorContBase::
VtkGridGeneratorContBaseSameDim<dim>::
VtkGridGeneratorContBaseSameDimCodim<codim>::
get_generators() -> GenMap_<dim, codim> &
{
  return grid_generators_;
}



Size
VtkIgaGridGeneratorContBase::
get_number_grids() const
{
  return generators_numbering_.size();
}



Size
VtkIgaGridGeneratorContBase::
get_number_active_grids() const
{
  Size count = 0;
  for (const auto &it : generators_numbering_)
    if (std::get<2>(it))
      ++count;

  return count;
}



const string &
VtkIgaGridGeneratorContBase::
get_grid_name(const Index &id) const
{
  Assert(id >= 0 && id < generators_numbering_.size(),
         ExcIndexRange(id, 0, generators_numbering_.size()));

  return std::get<1>(generators_numbering_[id]);
}



bool
VtkIgaGridGeneratorContBase::
get_grid_status(const std::string &name) const
{
  for (const auto &it : generators_numbering_)
  {
    if (std::get<1>(it) == name)
      return std::get<2>(it);
  }

  Assert(false, ExcMessage("Name not present."));
  return false; // Just for avoiding the compiler warning.
}



void
VtkIgaGridGeneratorContBase::
set_grid_status(const std::string &name, const bool status)
{
  for (auto &it : generators_numbering_)
  {
    if (std::get<1>(it) == name)
    {
      std::get<2>(it) = status;
      return;
    }
  }

  Assert(false, ExcMessage("Name not present."));
}



void
VtkIgaGridGeneratorContBase::
set_solid_grids(vtkMultiBlockDataSet *const mb)
{
  Assert(this->get_number_active_grids() > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(this->get_number_active_grids());

  unsigned int block_index = 0;
  boost::fusion::for_each(data_varying_dim_,
                          [&](const auto & type_and_data_same_dim)
  {
    boost::fusion::for_each(type_and_data_same_dim.second.get_data(),
                            [&](const auto & type_and_data_same_dim_codim)
    {
      auto &generators = type_and_data_same_dim_codim.second.get_generators();

      for (auto &g : generators)
      {
        const auto &grid_id = g.first;

        Assert(grid_id >= 0 && grid_id < generators_numbering_.size(),
               ExcIndexRange(grid_id, 0, generators_numbering_.size()));
        const auto &g_num = generators_numbering_[grid_id];

        // Is the grid active?
        if (std::get<2>(g_num))
        {
          const auto &name = std::get<1>(g_num);
          mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
          mb->SetBlock(block_index, g.second->get_solid_grid());
          ++block_index;
        }
      }

    } // end lambda function on codim
                           ); // for_each (data_varying_codim_

  } // end lambda function on dim
                         ); // for_each (data_varying_dim_
}



void
VtkIgaGridGeneratorContBase::
set_knot_grids(vtkMultiBlockDataSet *const mb)
{
  Assert(this->get_number_active_grids() > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(this->get_number_active_grids());

  Index block_index = 0;
  boost::fusion::for_each(data_varying_dim_,
                          [&](const auto & type_and_data_same_dim)
  {
    boost::fusion::for_each(type_and_data_same_dim.second.get_data(),
                            [&](const auto & type_and_data_same_dim_codim)
    {
      auto &generators = type_and_data_same_dim_codim.second.get_generators();

      for (auto &g : generators)
      {
        const auto &grid_id = g.first;

        Assert(grid_id >= 0 && grid_id < generators_numbering_.size(),
               ExcIndexRange(grid_id, 0, generators_numbering_.size()));
        const auto &g_num = generators_numbering_[grid_id];

        // Is the grid active?
        if (std::get<2>(g_num))
        {
          const auto &name = std::get<1>(g_num);
          mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
          mb->SetBlock(block_index, g.second->get_knot_grid());
          ++block_index;
        }
      }

    } // end lambda function on codim
                           ); // for_each (data_varying_codim_

  } // end lambda function on dim
                         ); // for_each (data_varying_dim_
}



VtkIgaGridGeneratorContParm::
VtkIgaGridGeneratorContParm(const FunContPtr_ funcs_container,
                            const GridInfoPtr_ solid_info,
                            const GridInfoPtr_ knot_info)
  :
  Base_(funcs_container, solid_info, knot_info, nullptr)
{
  this->fill_generators();
}



auto
VtkIgaGridGeneratorContParm::
create(const FunContPtr_ funcs_container,
       const GridInfoPtr_ solid_info,
       const GridInfoPtr_ knot_info) -> SelfPtr_
{
  return SelfPtr_(new Self_(funcs_container, solid_info, knot_info));
}


#if 0
void
VtkIgaGridGeneratorContParm::
update_parametric(const GridInfoPtr_ solid_info,
                  const GridInfoPtr_ knot_info)
{
  const bool solid_updated = solid_info_->update(solid_info);
  const bool knot_updated = knot_info_->update(knot_info);


  boost::fusion::for_each(data_varying_dim_,
                          [&](const auto & type_and_data_same_dim)
  {

    boost::fusion::for_each(type_and_data_same_dim.second.get_data(),
                            [&](const auto & type_and_data_same_dim_codim)
    {
      AssertThrow(false,ExcNotImplemented());
      //TODO: (martinelli,23 Oct 2015): the next commented loop causes internal compiler error on gcc-5.2.0
      /*
            auto &generators = type_and_data_same_dim_codim.second.get_generators();
            for (auto &g : generators)
              g.second->update(solid_updated, knot_updated,false);
            //*/
    } // end lambda function on codim
                           ); // for_each (data_varying_codim_

  } // end lambda function on dim
                         ); // for_each (data_varying_dim_
}
#endif

#if 0
void
VtkIgaGridGeneratorContParm::
fill_generators()
{
  // Iterating over all the functions in the container for building the generators
  // and setting up the numbering for them.

  const auto &funcs_container_data = funcs_container_->get_data();
  boost::fusion::for_each(
    funcs_container_data,
    [&](const auto & type_and_data_same_dim)
  {

    boost::fusion::for_each(type_and_data_same_dim.second.get_data(),
                            [&](const auto & type_and_data_same_dim_codim)
    {

      const auto &domains = type_and_data_same_dim_codim.second.get_all_domains();

      for (const auto &domain : domains)
      {
        this->insert_generator(domain.get_ptr_const_data());
        /*
        // Is it an identity mapping?
        if (std::dynamic_pointer_cast<IdFun>(map_fun) != nullptr)
          this->insert_generator<dim, codim> (map_fun, name);
        //*/
      } // endl loop on map_funs with a given pair <dim,codim>

    } // end lambda function on codim
                           ); // for_each (data_varying_codim_

  } // end lambda function on dim
  ); // for_each (data_varying_dim_
}
#endif

#if 0
template<int dim, int codim>
void
VtkIgaGridGeneratorContParm::
insert_generator(const DomainPtr_<dim, codim> domain)
{
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
                       (domain, solid_info_, knot_info_, funcs_container_);
}
#endif


VtkIgaGridGeneratorContPhys::
VtkIgaGridGeneratorContPhys(const FunContPtr_ funcs_container,
                            const GridInfoPtr_ solid_info,
                            const GridInfoPtr_ knot_info,
                            const ControlGridInfoPtr_ control_info)
  :
  Base_(funcs_container, solid_info, knot_info,control_info)
{
  Assert(control_info != nullptr, ExcNullPtr())
  this->fill_generators();
}



auto
VtkIgaGridGeneratorContPhys::
create(const FunContPtr_ funcs_container,
       const GridInfoPtr_ solid_info,
       const GridInfoPtr_ knot_info,
       const ControlGridInfoPtr_ control_info) -> SelfPtr_
{
  return SelfPtr_(new Self_(funcs_container, solid_info, knot_info, control_info));
}



void
VtkIgaGridGeneratorContBase::
update(const GridInfoPtr_ solid_info,
       const GridInfoPtr_ knot_info,
       const ControlGridInfoPtr_ control_info)
{
  const bool solid_updated = solid_info_->update(solid_info);
  const bool knot_updated = knot_info_->update(knot_info);

  bool control_updated;
  if (control_info_)
    control_updated = control_info_->update(control_info);
  else
    false;

  boost::fusion::for_each(data_varying_dim_,
                          [&](const auto & type_and_data_same_dim)
  {

    boost::fusion::for_each(type_and_data_same_dim.second.get_data(),
                            [&](const auto & type_and_data_same_dim_codim)
    {
      AssertThrow(false,ExcNotImplemented());

      //TODO: (martinelli,23 Oct 2015): the next commented loop causes internal compiler error on gcc-5.2.0

      auto &generators = type_and_data_same_dim_codim.second.get_generators();
      for (auto &g : generators)
        g.second->update(solid_updated, knot_updated, control_updated);
      //*/
    } // end lambda function on codim
                           ); // for_each (data_varying_codim_

  } // end lambda function on dim
                         ); // for_each (data_varying_dim_
}



void
VtkIgaGridGeneratorContPhys::
set_control_grids(vtkMultiBlockDataSet *const mb)
{
  Assert(this->get_number_active_grids_ig() > 0, ExcEmptyObject());

  mb->SetNumberOfBlocks(this->get_number_active_grids_ig());

  Index block_index = 0;
  boost::fusion::for_each(data_varying_dim_,
                          [&](const auto & type_and_data_same_dim)
  {
    using Type_Value = typename std::remove_reference<decltype(type_and_data_same_dim)>::type;
    using Type = typename Type_Value::first_type;
    const int dim = Type::value;

    boost::fusion::for_each(type_and_data_same_dim.second.get_data(),
                            [&](const auto & type_and_data_same_dim_codim)
    {
      using Type_Value = typename std::remove_reference<decltype(type_and_data_same_dim_codim)>::type;
      using Type = typename Type_Value::first_type;
      const int codim = Type::value;

      auto &generators = type_and_data_same_dim_codim.second.get_generators();

      for (auto &g : generators)
      {
        const auto &grid_id = g.first;

        Assert(grid_id >= 0 && grid_id < generators_numbering_.size(),
               ExcIndexRange(grid_id, 0, generators_numbering_.size()));
        const auto &g_num = generators_numbering_[grid_id];

        // Is the grid active?
        if (std::get<2>(g_num))
        {
          // Is an ig mapping?
          if (std::get<3>(g_num))
          {
            const auto &name = std::get<1>(g_num);
            mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), name.c_str());
            mb->SetBlock(block_index,
                         std::dynamic_pointer_cast<VtkIgaGridGeneratorPhys<dim,codim> >(
                           g.second)->get_control_grid());
            ++block_index;
          }

        }
      }

    } // end lambda function on codim
                           ); // for_each (data_varying_codim_

  } // end lambda function on dim
                         ); // for_each (data_varying_dim_
}



Size
VtkIgaGridGeneratorContPhys::
get_number_active_grids_ig() const
{
  Size count = 0;
  for (const auto &it : generators_numbering_)
    if (std::get<3>(it))
      ++count;

  return count;
}



void
VtkIgaGridGeneratorContBase::
fill_generators()
{
  // Iterating over all the functions in the container for building the generators
  // and setting up the numbering for them.

  const auto &funcs_container_data = funcs_container_->get_data();
  boost::fusion::for_each(
    funcs_container_data,
    [&](const auto & type_and_data_same_dim)
  {
    boost::fusion::for_each(type_and_data_same_dim.second.get_data(),
                            [&](const auto & type_and_data_same_dim_codim)
    {
      const auto &domains = type_and_data_same_dim_codim.second.get_all_domains();

      for (const auto &domain : domains)
      {
        this->insert_generator(domain.get_ptr_const_data());
        // Is it a physical mapping?
        /*
        if (std::dynamic_pointer_cast<IdFun>(map_fun) == nullptr)
          this->insert_generator<dim, codim> (map_fun, name);
        //*/
      } // endl loop on map_funs with a given pair <dim,codim>

    } // end lambda function on codim
                           ); // for_each (data_varying_codim_

  } // end lambda function on dim
  ); // for_each (data_varying_dim_
}



template<int dim, int codim>
void
VtkIgaGridGeneratorContBase::
insert_generator(const DomainPtr_<dim, codim> domain)
{
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
                         (domain, solid_info_, knot_info_, control_info_, funcs_container_);

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
                         (domain, solid_info_, knot_info_, funcs_container_);
  }
}

IGA_NAMESPACE_CLOSE
