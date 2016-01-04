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

#include <paraview_plugin/grid_generator.h>

#include <igatools/geometry/domain.h>

#include <paraview_plugin/grid_information.h>
#include <paraview_plugin/knot_grid_generator.h>
#include <paraview_plugin/solid_grid_generator.h>
#include <paraview_plugin/control_grid_generator.h>


IGA_NAMESPACE_OPEN

template <class Domain>
VtkIgaGridGenerator<Domain>::
VtkIgaGridGenerator(const DomainPtr_ domain,
                    const Index &id,
                    const GridInfoPtr_ solid_grid_info,
                    const GridInfoPtr_ knot_grid_info,
                    const ControlGridInfoPtr_ control_grid_info,
                    const ObjContPtr_ obj_container,
                    const bool is_active,
                    const bool is_physical)
  :
  domain_(domain),
  id_ (id),
  solid_grid_info_(solid_grid_info),
  knot_grid_info_(knot_grid_info),
  control_grid_info_(control_grid_info),
  objs_container_(obj_container),
  solid_grid_(VtkGridPtr_()),
  knot_grid_(VtkGridPtr_()),
  control_grid_(VtkGridPtr_()),
  recompute_solid_(true),
  recompute_knot_(true),
  recompute_control_(is_physical ? true : false),
  is_active_ (is_active),
  is_physical_ (is_physical),
  is_ig_grid_func_ (is_physical ? std::dynamic_pointer_cast<const IgGridFunc_>(domain->get_grid_function()) != nullptr  : false)
{
  Assert(domain_ != nullptr, ExcNullPtr());
  Assert(solid_grid_info_ != nullptr, ExcNullPtr());
  Assert(knot_grid_info_ != nullptr, ExcNullPtr());
  Assert(control_grid_info_ != nullptr, ExcNullPtr());
  Assert(objs_container_ != nullptr, ExcNullPtr());
}



template <class Domain>
auto
VtkIgaGridGenerator<Domain>::
create(const DomainPtr_ domain,
       const Index &id,
       const GridInfoPtr_ solid_grid_info,
       const GridInfoPtr_ knot_grid_info,
       const ControlGridInfoPtr_ control_grid_info,
       const ObjContPtr_ obj_container,
       const bool is_active,
       const bool is_physical) -> SelfPtr_
{
    return SelfPtr_ (new Self_(domain, id, solid_grid_info,
                               knot_grid_info, control_grid_info,
                               obj_container, is_active, is_physical));
}



template <class Domain>
const Index &
VtkIgaGridGenerator<Domain>::
get_id() const
{
    return id_;
}



template <class Domain>
const std::string &
VtkIgaGridGenerator<Domain>::
get_name() const
{
    return domain_->get_name();
}



template <class Domain>
bool
VtkIgaGridGenerator<Domain>::
is_active() const
{
    return is_active_;
}



template <class Domain>
void
VtkIgaGridGenerator<Domain>::
set_status(const bool status_flag)
{
    is_active_ = status_flag;
}



template <class Domain>
bool
VtkIgaGridGenerator<Domain>::
is_ig_grid_func() const
{
    return is_ig_grid_func_;
}



template <class Domain>
bool
VtkIgaGridGenerator<Domain>::
is_physical() const
{
    return is_physical_;
}



template <class Domain>
auto
VtkIgaGridGenerator<Domain>::
get_solid_grid() -> VtkGridPtr_
{
  if (recompute_solid_)
  {
    // Recomputing solid grid.
    solid_grid_ = VtkIgaSolidGridGenerator<Domain>::
            create_grid(domain_, solid_grid_info_, objs_container_);

    recompute_solid_ = false;
  }

  Assert(solid_grid_ != nullptr, ExcNullPtr());
  return solid_grid_;
}



template <class Domain>
auto
VtkIgaGridGenerator<Domain>::
get_knot_grid() -> VtkGridPtr_
{
  if (recompute_knot_)
  {
    // Recomputing knot grid.
    knot_grid_ = VtkIgaKnotGridGenerator<Domain>::
            create_grid(domain_, knot_grid_info_);

    recompute_knot_ = false;
  }

  Assert(knot_grid_ != nullptr, ExcNullPtr());
  return knot_grid_;
}



template <class Domain>
auto
VtkIgaGridGenerator<Domain>::
get_control_grid() -> VtkGridPtr_
{
  Assert (is_physical_,
          ExcMessage("Control mesh cannot be retrieved for a parametric "
                     "mesh."));

  if (recompute_control_)
  {
    // Recomputing control grid.
    control_grid_ = VtkIgaControlGridGenerator<Domain>::
            create_grid(this->domain_, control_grid_info_);

    recompute_control_ = false;
  }

  Assert(control_grid_ != nullptr, ExcNullPtr());
  return control_grid_;
}



template <class Domain>
void
VtkIgaGridGenerator<Domain>::
update(const bool solid_updated,
       const bool knot_updated,
       const bool control_updated)
{
  if (!this->recompute_solid_)
    this->recompute_solid_ = solid_updated;

  if (!this->recompute_knot_)
    this->recompute_knot_ = knot_updated;

  Assert (is_physical_ || !control_updated,
          ExcMessage("Control mesh cannot be updated for a parametric "
                     "mesh."));
  if (!this->recompute_control_ && is_physical_)
    this->recompute_control_ = control_updated;
}


template class VtkIgaGridGenerator<Domain<1, 0>>;
template class VtkIgaGridGenerator<Domain<1, 1>>;
template class VtkIgaGridGenerator<Domain<1, 2>>;
template class VtkIgaGridGenerator<Domain<2, 0>>;
template class VtkIgaGridGenerator<Domain<2, 1>>;
template class VtkIgaGridGenerator<Domain<3, 0>>;


IGA_NAMESPACE_CLOSE
