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

#include <igatools/geometry/identity_mapping.h>
#include <igatools/base/exceptions.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>

using std::vector;
using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

template<int dim, int codim>
IdentityMapping< dim, codim>::
IdentityMapping(const std::shared_ptr<GridType> grid)
    :
    Mapping<dim, codim>(grid)
{
    for (int i = 0; i < dim; ++i)
        A_[i][i] = 1.;

    for (Index face_id = 0; face_id < UnitElement<dim>::faces_per_element; ++face_id)
    {
        Index k = 0 ;
        for (auto &i : UnitElement<dim>::face_active_directions[face_id])
            face_A_[face_id][k++][i] = 1. ;
    }

}



template<int dim, int codim>
IdentityMapping< dim, codim>::
IdentityMapping(const IdentityMapping<dim,codim> &map)
    :
    Mapping<dim,codim>::Mapping(map)
{}



template<int dim, int codim>
auto
IdentityMapping< dim, codim>::
create(const std::shared_ptr<GridType> grid) -> shared_ptr<base_t>
{
    return (shared_ptr<base_t>(new self_t(grid)));
}



template<int dim, int codim>
void
IdentityMapping< dim, codim>::
evaluate(vector<ValueType> &values) const
{
    const int num_points = points_.size();
    Assert(values.size() == num_points,
           ExcDimensionMismatch(values.size(),num_points));
    for (int i = 0; i<num_points; i++)
    {
        const auto &x = points_[i];
        for (int k = 0; k < dim; ++k)
        {
            values[i][k] = x[k];
        }
        for (int k = dim; k < codim; ++k)
        {
            values[i][k] = 0.;
        }
    }
}



template<int dim, int codim>
void
IdentityMapping< dim, codim>::
evaluate_gradients(vector<GradientType> &gradients) const
{
    Assert(gradients.size() == points_.size(),
           ExcDimensionMismatch(gradients.size(),points_.size()));

    const int num_points = points_.size();

    for (int i = 0; i<num_points; i++)
        gradients[i] = A_;
}



template<int dim, int codim>
void
IdentityMapping< dim, codim>::
evaluate_hessians(vector<HessianType> &hessians) const
{
    const int num_points = points_.size();
    for (int i = 0; i<num_points; i++)
        hessians[i] = 0.;
}



template<int dim, int codim>
void
IdentityMapping< dim, codim>::
evaluate_face(const Index face_id, vector<ValueType> &values) const
{
    Assert(face_id < UnitElement<dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim>::faces_per_element));
    const auto &points = face_points_[face_id] ;
    const int num_points = points.size();
    Assert(values.size() == num_points,
           ExcDimensionMismatch(values.size(),num_points));
    for (int i = 0; i<num_points; i++)
    {
        const auto &x = points[i];
        for (int k = 0; k < dim; ++k)
        {
            values[i][k] = x[k];
        }
        for (int k = dim; k < codim; ++k)
        {
            values[i][k] = 0.;
        }
    }
}



template<int dim, int codim>
void
IdentityMapping< dim, codim>::
evaluate_face_gradients(const Index face_id, vector<GradientFaceType> &gradients) const
{
    Assert(face_id < UnitElement<dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim>::faces_per_element));
    const int num_points = face_points_[face_id].size();
    Assert(gradients.size() == num_points,
           ExcDimensionMismatch(gradients.size(),num_points));


    for (int i = 0; i<num_points; i++)
        gradients[i] = face_A_[face_id];
}



template<int dim, int codim>
void
IdentityMapping< dim, codim>::
evaluate_face_hessians(const Index face_id, vector<HessianFaceType> &hessians) const
{
    Assert(face_id < UnitElement<dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim>::faces_per_element));
    const int num_points = face_points_[face_id].size();
    for (int i = 0; i<num_points; i++)
        hessians[i] = 0.;
}



template<int dim, int codim>
ValueFlags
IdentityMapping< dim, codim>::
required_flags() const
{
    return ValueFlags::point;
}



template<int dim, int codim>
void
IdentityMapping< dim, codim>::
set_element(const CartesianGridElementAccessor<dim> &elem)
{
    points_ = elem.get_points();
}



template<int dim, int codim>
void
IdentityMapping< dim, codim>::
set_face_element(const Index face_id,
                 const CartesianGridElementAccessor<dim> &elem)
{
    face_points_[face_id] = elem.get_face_points(face_id);
}



template<int dim, int codim>
shared_ptr< Mapping< dim, codim> >
IdentityMapping< dim, codim>::
clone() const
{
    return shared_ptr<self_t>(new self_t(*this));
}



template<int dim, int codim>
void
IdentityMapping< dim, codim>::
print_info(LogStream &out) const
{
    out << "Type = IdentityMapping<" << dim;
    out << "," << codim << ">" << std::endl;
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/identity_mapping.inst>
