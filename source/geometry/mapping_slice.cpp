//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

#include <igatools/geometry/mapping_slice.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/exceptions.h>

using std::array;
using std::vector;
using std::shared_ptr;
using std::endl;

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
MappingSlice<dim_, codim_>::
MappingSlice(const std::shared_ptr<const SupMap> map,
             const int face_id,
             const std::shared_ptr<GridType > grid,
             const std::shared_ptr<std::map<int,int> > elem_map)
    :
    base_t::Mapping(grid),
    map_(map),
    direction_(UnitElement<dim + 1>::face_constant_direction[face_id]),
    value_(UnitElement<dim + 1>::face_side[face_id]),
    element(map_->begin()),
    elem_map_ {elem_map}
{}



template<int dim_, int codim_>
MappingSlice<dim_, codim_>::
MappingSlice(const self_t &map_slice)
    :
    base_t::Mapping(map_slice),
    map_(map_slice.map_),
    direction_(map_slice.direction_),
    value_(map_slice.value_),
    element(map_slice.element)
{}



template<int dim_, int codim_>
auto
MappingSlice<dim_, codim_>::
create(const std::shared_ptr<const SupMap> map,
       const int face_id,
       const std::shared_ptr<GridType > grid,
       const std::shared_ptr<std::map<int,int> > elem_map) -> shared_ptr<base_t>
{
    return shared_ptr<base_t>(new self_t(map, face_id, grid, elem_map));
}



template<int dim_, int codim_>
auto
MappingSlice<dim_, codim_>::
build_extended_quadrature(const Quadrature<dim> &quad) const -> Quadrature<dim+1>
{
    const auto points  = quad.get_points();
    const auto weights = quad.get_weights();

    auto ext_quad = Quadrature<dim+1>(
                        insert(points, direction_,std::vector<Real>(1,value_)),
                        insert(weights,direction_,std::vector<Real>(1,1.0))) ;

    return ext_quad;
}



template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
evaluate(std::vector<Value> &values) const
{
    values = element->get_map_values();
}



template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
evaluate_gradients(std::vector<Gradient> &gradients) const
{
    auto grad = element->get_map_gradients();

    const auto active_dir = UnitElement<dim+1>::active_directions[direction_];

    const int num_points = grad.size();
    for (int p = 0 ; p < num_points; p++)
        for (int i = 0; i < dim; i++)
            gradients[p][i] = grad[p][active_dir[i]];

}



template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
init_element(const ValueFlags flag, const Quadrature<dim> &quad) const
{
    element->init_cache(flag, build_extended_quadrature(quad));
}



template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
set_element(const CartesianGridElementAccessor<dim> &elem) const
{
    element->move_to((*elem_map_)[elem.get_flat_index()]);
    element->fill_cache();
}


template<int dim_, int codim_>
void
MappingSlice<dim_,codim_>::
set_face_element(const Index face_id,
                 const CartesianGridElementAccessor<dim> &elem) const
{
    Assert(false, ExcNotImplemented());
    AssertThrow(false, ExcNotImplemented());
}


//TODO(pauletti, Jun 20, 2014): simplify the output
template<int dim_, int codim_>
void
MappingSlice<dim_, codim_>::
print_info(LogStream &out) const
{
    out << "Type = MappingSlice<" << dim_ << "," << dim_+codim_ << ">" << endl;

    out.push("\t");
    out << "Direction = " << direction_ << endl ;
    out << "    ValueType = " << value_ << endl ;

    out << "Sliced Map:" << endl ;
    out.push("\t");
    map_->print_info(out) ;
    out << endl;

    out.pop();
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/mapping_slice.inst>
