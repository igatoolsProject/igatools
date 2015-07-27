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

#include <paraview_plugin/control_grid_generator.h>
#include <paraview_plugin/grid_information.h>
#include <paraview_plugin/knot_grid_generator.h>
#include <paraview_plugin/solid_grid_generator.h>


IGA_NAMESPACE_OPEN

template <int dim, int codim>
VtkIgaGridGenerator<dim, codim>::
VtkIgaGridGenerator (const MapFunPtr_ map_fun,
                  const GridInfoPtr_ solid_grid_info,
                  const GridInfoPtr_ knot_grid_info,
                  const ControlGridInfoPtr_ control_grid_info,
                  const FunContPtr_ func_container,
                  const bool is_physical)
    :
                  map_fun_ (map_fun),
                  solid_grid_info_ (solid_grid_info),
                  knot_grid_info_ (knot_grid_info),
                  control_grid_info_ (control_grid_info),
                  funcs_container_ (func_container),
                  is_physical_ (is_physical),
                  recompute_solid_ (true),
                  recompute_knot_ (true),
                  recompute_control_ (true)
{
    Assert (map_fun != nullptr, ExcNullPtr());
    Assert (solid_grid_info_ != nullptr, ExcNullPtr());
    Assert (knot_grid_info_ != nullptr, ExcNullPtr());

    Assert (func_container != nullptr, ExcNullPtr());
}



template <int dim, int codim>
auto
VtkIgaGridGenerator<dim, codim>::
create_physical (const MapFunPtr_ map_fun,
                 const GridInfoPtr_ solid_grid_info,
                 const GridInfoPtr_ knot_grid_info,
                 const ControlGridInfoPtr_ control_grid_info,
                 const FunContPtr_ func_container) -> SelfPtr_
{
    Assert (control_grid_info != nullptr, ExcNullPtr());

    return SelfPtr_ (new Self_(map_fun, solid_grid_info, knot_grid_info,
                               control_grid_info, func_container, true));
}



template <int dim, int codim>
auto
VtkIgaGridGenerator<dim, codim>::
create_parametric (const MapFunPtr_ map_fun,
                   const GridInfoPtr_ solid_grid_info,
                   const GridInfoPtr_ knot_grid_info,
                   const FunContPtr_ func_container) -> SelfPtr_
{
    return SelfPtr_ (new Self_(map_fun, solid_grid_info, knot_grid_info,
                               ControlGridInfoPtr_(), func_container, false));
}



template <int dim, int codim>
void
VtkIgaGridGenerator<dim, codim>::
update (const bool solid_updated,
        const bool knot_updated,
        const bool control_updated)
{
    Assert (is_physical_, ExcMessage("Not valid for parametric grids."));

    if (!recompute_solid_)
        recompute_solid_ = solid_updated;

    if (!recompute_knot_)
        recompute_knot_ = knot_updated;

    if (!recompute_control_)
        recompute_control_ = control_updated;

}

template <int dim, int codim>
void
VtkIgaGridGenerator<dim, codim>::
update (const bool solid_updated,
        const bool knot_updated)
{
    Assert (!is_physical_, ExcMessage("Not valid for physical grids."));

    if (!recompute_solid_)
        recompute_solid_ = solid_updated;

    if (!recompute_knot_)
        recompute_knot_ = knot_updated;
}



template <int dim, int codim>
auto
VtkIgaGridGenerator<dim, codim>::
get_solid_grid () -> VtkGridPtr_
{
    if (recompute_solid_)
    {
        solid_grid_ = VtkIgaSolidGridGenerator<dim, codim>::
                get_grid(map_fun_, solid_grid_info_, funcs_container_);

        recompute_solid_ = false;
    }

    Assert (solid_grid_ != nullptr, ExcNullPtr());
    return solid_grid_;
}



template <int dim, int codim>
auto
VtkIgaGridGenerator<dim, codim>::
get_knot_grid () -> VtkGridPtr_
{
    if (recompute_knot_)
    {
        knot_grid_ = VtkIgaKnotGridGenerator<dim, codim>::
                get_grid(map_fun_, knot_grid_info_);

        recompute_knot_ = false;
    }

    Assert (knot_grid_ != nullptr, ExcNullPtr());
    return knot_grid_;
}



template <int dim, int codim>
auto
VtkIgaGridGenerator<dim, codim>::
get_control_grid () -> VtkGridPtr_
{
    Assert (is_physical_, ExcMessage("Not valid for parametric grids."));

    if (recompute_control_)
    {
        control_grid_ = VtkIgaControlGridGenerator<dim, codim>::
                get_grid(map_fun_, control_grid_info_);

        recompute_control_ = false;
    }

    Assert (control_grid_ != nullptr, ExcNullPtr());
    return control_grid_;
}



template class VtkIgaGridGenerator<1, 0>;
template class VtkIgaGridGenerator<1, 1>;
template class VtkIgaGridGenerator<1, 2>;
template class VtkIgaGridGenerator<2, 0>;
template class VtkIgaGridGenerator<2, 1>;
template class VtkIgaGridGenerator<3, 0>;



IGA_NAMESPACE_CLOSE
