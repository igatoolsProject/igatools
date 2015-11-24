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

#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/NURBS.h>

#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/container/map/map_fwd.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <boost/fusion/include/at_key.hpp>


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
insert_grid (const GridPtr_t<dim> grid, const Index &id)
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
insert_spline_space (const SplineSpacePtr_t<dim, range, rank> spline_space,
                  const Index &id)
{
    Assert (spline_space != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert ((!this->is_spline_space<dim, range, rank>(id)),
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
    Assert ((this->is_spline_space<dim, range, rank>(id)),
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
insert_ref_space (const RefSpacePtr_t<dim, range, rank> ref_space,
               const Index &id)
{
    Assert (ref_space != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert ((!this->is_ref_space<dim, range, rank>(id)),
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
    Assert ((this->is_ref_space<dim, range, rank>(id)),
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



template <int dim, int range, int rank>
auto
ObjectsContainer::
get_bspline (const Index &id) const -> shared_ptr<BSpline<dim, range, rank>>
{
    Assert ((this->template is_bspline<dim, range, rank>(id)),
            ExcMessage("Object id does not correspond to a BSpline "
                       "space with the given dimensions"));
    const auto ref_space = at_key<RefSpace_t<dim, range, rank>>(objects_).at(id);
    auto bspline = std::dynamic_pointer_cast<BSpline<dim, range, rank>> (ref_space);
    Assert (bspline != nullptr, ExcNullPtr());
    return bspline;
};



template <int dim, int range, int rank>
auto
ObjectsContainer::
get_nurbs (const Index &id) const -> shared_ptr<NURBS<dim, range, rank>>
{
    Assert ((this->template is_nurbs<dim, range, rank>(id)),
            ExcMessage("Object id does not correspond to a NURBS "
                       "space with the given dimensions"));
    const auto ref_space = at_key<RefSpace_t<dim, range, rank>>(objects_).at(id);
    auto nurbs = std::dynamic_pointer_cast<NURBS<dim, range, rank>> (ref_space);
    Assert (nurbs != nullptr, ExcNullPtr());
    return nurbs;
};



template <int dim, int range, int rank>
bool
ObjectsContainer::
is_bspline (const Index &id) const
{
    if (this->is_ref_space<dim, range, rank>(id))
    {
        const auto ref_space = at_key<RefSpace_t<dim, range, rank>>(objects_).at(id);
        if (std::dynamic_pointer_cast<BSpline<dim, range, rank>>(ref_space) != nullptr)
            return true;
        else
            return false;
    }
    else
        return false;
};



template <int dim, int range, int rank>
bool
ObjectsContainer::
is_nurbs (const Index &id) const
{
    if (this->is_ref_space<dim, range, rank>(id))
    {
        const auto ref_space = at_key<RefSpace_t<dim, range, rank>>(objects_).at(id);
        if (std::dynamic_pointer_cast<NURBS<dim, range, rank>>(ref_space) != nullptr)
            return true;
        else
            return false;
    }
    else
        return false;
};



template <int dim, int range>
void
ObjectsContainer::
insert_grid_function (const GridFuncPtr_t<dim, range> grid_func,
                   const Index &id)
{
    Assert (grid_func != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert ((!this->is_grid_function<dim, range>(id)),
            ExcMessage("Object id already defined for a grid function "
                       "with the given dimensions"));

    auto &grid_funcs = at_key<GridFunction<dim, range>>(objects_);

    grid_funcs[id] = grid_func;
};



template <int dim, int range>
auto
ObjectsContainer::
get_grid_function (const Index &id) const -> GridFuncPtr_t<dim, range>
{
    Assert ((this->is_grid_function<dim, range>(id)),
            ExcMessage("Object id does not correspond to a grid function "
                       "with the given dimensions"));
    return at_key<GridFunction<dim, range>>(objects_).at(id);
};



template <int dim, int range>
bool
ObjectsContainer::
is_grid_function (const Index &id) const
{
    const auto &grid_funcs = at_key<GridFunction<dim, range>>(objects_);
    return grid_funcs.find(id) != grid_funcs.end();
};



template <int dim, int codim>
void
ObjectsContainer::
insert_domain (const DomainPtr_t<dim, codim> domain,
            const Index &id)
{
    Assert (domain != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert ((!this->is_domain<dim, codim>(id)),
            ExcMessage("Object id already defined for a domain "
                       "with the given dimensions"));

    auto &domains = at_key<Domain<dim, codim>>(objects_);

    domains[id] = domain;
};



template <int dim, int codim>
auto
ObjectsContainer::
get_domain (const Index &id) const -> DomainPtr_t<dim, codim>
{
    Assert ((this->is_domain<dim, codim>(id)),
            ExcMessage("Object id does not correspond to a domain "
                       "with the given dimensions"));
    return at_key<Domain<dim, codim>>(objects_).at(id);
};



template <int dim, int codim>
bool
ObjectsContainer::
is_domain (const Index &id) const
{
    const auto &domains = at_key<Domain<dim, codim>>(objects_);
    return domains.find(id) != domains.end();
};



template <int dim, int range, int rank, int codim>
void
ObjectsContainer::
insert_phys_space_basis (const PhysSpacePtr_t<dim, range, rank, codim> space,
                      const Index &id)
{
    Assert (space != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert ((!this->is_phys_space_basis<dim, range, rank, codim>(id)),
            ExcMessage("Object id already defined for a physical space "
                       "basis with the given dimensions"));

    auto &spaces = at_key<PhysSpace_t<dim, range, rank, codim>>(objects_);

    spaces[id] = space;
};



template <int dim, int range, int rank, int codim>
auto
ObjectsContainer::
get_phys_space_basis (const Index &id) const ->
PhysSpacePtr_t<dim, range, rank, codim>
{
    Assert ((this->is_phys_space_basis<dim, range, rank, codim>(id)),
            ExcMessage("Object id does not correspond to a physical space"
                       " basis with the given dimensions"));
    return at_key<PhysSpace_t<dim, range, rank, codim>>(objects_).at(id);
};



template <int dim, int range, int rank, int codim>
bool
ObjectsContainer::
is_phys_space_basis (const Index &id) const
{
    const auto &spaces = at_key<PhysSpace_t<dim, range, rank, codim>>(objects_);
    return spaces.find(id) != spaces.end();
};



template <int dim, int codim, int range, int rank>
void
ObjectsContainer::
insert_function (const FunctionPtr_t<dim, codim, range, rank> function,
              const Index &id)
{
    Assert (function != nullptr, ExcNullPtr());

    // TODO: this will be replaced when the global enumeration control is
    // performed.
    Assert ((!this->is_function<dim, codim, range, rank>(id)),
            ExcMessage("Object id already defined for a function "
                       "with the given dimensions"));

    auto &functions = at_key<Function<dim, codim, range, rank>>(objects_);

    functions[id] = function;
};



template <int dim, int codim, int range, int rank>
auto
ObjectsContainer::
get_function (const Index &id) const ->
FunctionPtr_t<dim, codim, range, rank>
{
    Assert ((this->is_function<dim, codim, range, rank>(id)),
            ExcMessage("Object id does not correspond to a function "
                       "with the given dimensions"));
    return at_key<Function<dim, codim, range, rank>>(objects_).at(id);
};



template <int dim, int codim, int range, int rank>
bool
ObjectsContainer::
is_function (const Index &id) const
{
    const auto &functions = at_key<Function<dim, codim, range, rank>>(objects_);
    return functions.find(id) != functions.end();
};



IGA_NAMESPACE_CLOSE

#include <igatools/base/objects_container.inst>
