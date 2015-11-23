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

#include <igatools/base/objects_container.h>

#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/container/map/map_fwd.hpp>
#include <boost/fusion/include/map_fwd.hpp>

using std::shared_ptr;
using namespace boost::fusion;

IGA_NAMESPACE_OPEN

ObjectsContainer::
ObjectsContainer()
{
};



auto
ObjectsContainer::
create() -> shared_ptr<self_t>
{
    return shared_ptr<self_t> (new self_t());
};



auto
ObjectsContainer::
const_create() -> shared_ptr<const self_t>
{
    return shared_ptr<const self_t> (new self_t());
};



template <int dim>
void
ObjectsContainer::
add_grid (const GridPtr_t<dim> grid, const Index &id)
{
    Assert (grid != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert (!this->is_grid<dim>(id),
            ExcMessage("Object id already defined for a grid with the "
                       "given dim"));

    auto &grids = at_key<Grid<dim>>(objects_);

    grids[id] = grid;
};



template <int dim>
auto
ObjectsContainer::
get_grid (const Index &id) const -> GridPtr_t<dim>
{
    Assert (this->is_grid<dim>(id),
            ExcMessage("Object id does not correspond to a grid with the "
                       "given dim"));
    return at_key<Grid<dim>>(objects_).at(id);
};



template <int dim>
bool
ObjectsContainer::
is_grid (const Index &id) const
{
    const auto &grids = at_key<Grid<dim>>(objects_);
    return grids.find(id) != grids.end();
};



template <int dim, int range, int rank>
void
ObjectsContainer::
add_spline_space (const SplineSpacePtr_t<dim, range, rank> spline_space,
                  const Index &id)
{
    Assert (spline_space != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert (!this->is_spline_space<dim>(id),
            ExcMessage("Object id already defined for a spline space with "
                       "the given dimensions"));

    auto &spline_spaces = at_key<SplineSpace<dim, range, rank>>(objects_);

    spline_spaces[id] = spline_space;
};



template <int dim, int range, int rank>
auto
ObjectsContainer::
get_spline_space (const Index &id) const -> SplineSpacePtr_t<dim, range, rank>
{
    Assert (this->is_spline_space<dim>(id),
            ExcMessage("Object id does not correspond to a spline space "
                       "with the given dimensions"));
    return at_key<SplineSpace<dim, range, rank>>(objects_).at(id);
};



template <int dim, int range, int rank>
bool
ObjectsContainer::
is_spline_space (const Index &id) const
{
    const auto &spline_spaces = at_key<SplineSpace<dim, range, rank>>(objects_);
    return spline_spaces.find(id) != spline_spaces.end();
};



template <int dim, int range, int rank>
void
ObjectsContainer::
add_ref_space (const RefSpacePtr_t<dim, range, rank> ref_space,
               const Index &id)
{
    Assert (ref_space != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert (!this->is_ref_space<dim>(id),
            ExcMessage("Object id already defined for a reference space "
                       "with the given dimensions"));

    auto &ref_spaces = at_key<RefSpace_t<dim, range, rank>>(objects_);

    ref_spaces[id] = ref_space;
};



template <int dim, int range, int rank>
auto
ObjectsContainer::
get_ref_space (const Index &id) const -> RefSpacePtr_t<dim, range, rank>
{
    Assert (this->is_ref_space<dim>(id),
            ExcMessage("Object id does not correspond to a reference "
                       "space with the given dimensions"));
    return at_key<RefSpace_t<dim, range, rank>>(objects_).at(id);
};



template <int dim, int range, int rank>
bool
ObjectsContainer::
is_ref_space (const Index &id) const
{
    const auto &ref_spaces = at_key<RefSpace_t<dim, range, rank>>(objects_);
    return ref_spaces.find(id) != ref_spaces.end();
};



IGA_NAMESPACE_CLOSE

// #include <igatools/base/objects_container.inst>

