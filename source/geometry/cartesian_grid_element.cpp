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


#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/geometry/unit_element.h>
#include <igatools/utils/multi_array_utils.h>

#include <algorithm>

using std::array;
using std::vector;

IGA_NAMESPACE_OPEN

//TODO: inline the appropriate method and put in separate file

template <int dim>
CartesianGridElement<dim>::
CartesianGridElement(
    const CartesianGrid<dim> &grid,
    const Index index)
    :
    grid_(&grid)
{
    this->reset_flat_tensor_indices(index);
}


template <int dim>
const CartesianGrid<dim> *
CartesianGridElement<dim>::
get_grid() const
{
    return grid_;
}


template <int dim>
inline
auto
CartesianGridElement<dim>::
get_flat_index() const -> Index
{
    return flat_index_ ;
}


template <int dim>
inline
auto
CartesianGridElement<dim>::
get_tensor_index() const -> TensorIndex<dim>
{
    return tensor_index_ ;
}



template <int dim>
void
CartesianGridElement<dim>::
reset_flat_tensor_indices(const Index flat_index)
{
    Assert((flat_index == IteratorState::pass_the_end) ||
           ((flat_index >= 0) && (flat_index < grid_->get_num_elements())),
           ExcIndexRange(flat_index, 0, grid_->get_num_elements()));

    flat_index_ = flat_index ;

    //Fill tensor_index_
    if (flat_index_ != IteratorState::pass_the_end)
    {
        using Utils = MultiArrayUtils<dim>;
        tensor_index_ = Utils::flat_to_tensor_index(
                            flat_index_,
                            Utils::compute_weight(grid_->get_num_elements_dim()));
    }
    else
        tensor_index_.fill(IteratorState::pass_the_end);
}

template <int dim>
void
CartesianGridElement<dim>::
reset_flat_tensor_indices(const TensorIndex<dim> &tensor_index)
{
    tensor_index_= tensor_index;

    using Utils = MultiArrayUtils<dim>;
    flat_index_ = Utils::tensor_to_flat_index(
                      tensor_index_,
                      Utils::compute_weight(grid_->get_num_elements_dim()));

    Assert((flat_index_ == IteratorState::pass_the_end) ||
           ((flat_index_ >= 0) && (flat_index_ < grid_->get_num_elements())),
           ExcIndexRange(flat_index_, 0, grid_->get_num_elements()));
}



template <int dim>
Point<dim>
CartesianGridElement<dim>::
vertex(const int i) const
{
    Assert(i < UnitElement<dim>::vertices_per_element,
           ExcIndexRange(i,0,UnitElement<dim>::vertices_per_element));

    TensorIndex<dim> index = this->get_tensor_index();

    for (int j = 0; j < dim; ++j)
        index[j] += UnitElement<dim>::vertex_to_component[i][j];


    return this->get_grid()->get_knot_coordinates().cartesian_product(index);
};


template <int dim>
array< Real, dim>
CartesianGridElement<dim>::
get_coordinate_lengths() const
{
    const Point<dim> &p_origin = this->vertex(0);
    const Point<dim> &p_end = this->vertex(UnitElement<dim>::vertices_per_element-1);

    array<Real,dim> coord_length;
    for (int d = 0; d < dim ; ++d)
        coord_length[d] = p_end[d] - p_origin[d];

    return coord_length;
}


template <int dim>
Point< dim >
CartesianGridElement<dim>::
center() const
{
    Point<dim> center(vertex(0));
    center += vertex(UnitElement<dim>::opposite_vertex[0]);
    center *= 0.5;

    return (center) ;
}


template <int dim>
Real
CartesianGridElement<dim>::
get_measure(const TopologyId<dim> &topology_id) const
{
	const auto active_directions = topology_id.get_active_directions();

    Real result = 1.;
    for (const int &d : active_directions)
        result *= this->get_coordinate_lengths()[d];

    return result;
}



template <int dim>
Real
CartesianGridElement<dim>::
get_face_measure(const Index face_id) const
{
    return this->get_measure(FaceTopology<dim>(face_id));
}


template <int dim>
bool
CartesianGridElement<dim>::
is_point_inside(const Point< dim > &point) const
{
    const auto &knots_directions = this->get_grid()->get_knot_coordinates();
    const auto &tensor_index = this->get_tensor_index();
    for (int j = 0; j < dim; ++j)
    {
        const auto &knots = knots_directions.get_data_direction(j) ;
        const Index id = tensor_index[j] ;

        const Real pt_coord = point(j) ;
        if (pt_coord <= knots[id] || pt_coord >= knots[id+1])
            return (false) ;
    }

    return (true) ;
}



template <int dim>
bool CartesianGridElement<dim>::
is_boundary() const
{
    const auto num_elements_dim = this->get_grid()->get_num_elements_dim();

    const auto &element_index = this->get_tensor_index() ;

    for (int i = 0; i < dim; ++i)
    {
        if (element_index[i] == 0 or element_index[i] == num_elements_dim(i) - 1)
            return true;
    }

    return false;
}


template <int dim>
bool
CartesianGridElement<dim>::
is_boundary(const Index face_id) const
{
    const int const_direction = UnitElement<dim>::face_constant_direction[face_id];
    const int face_side = UnitElement<dim>::face_side[face_id];

    const auto element_id_dir = this->get_tensor_index()[const_direction] ;
    const auto num_elements_dir = this->get_grid()->get_num_elements_dim()(const_direction);

    return (element_id_dir == ((num_elements_dir-1) * face_side)) ;
}


template <int dim>
void
CartesianGridElement<dim>::
print_info(LogStream &out,const VerbosityLevel verbosity) const
{
    using std::endl;

    const std::string tab = "   ";

    out << "CartesianGridElement<" << dim << "> info:" << endl;
    out.push(tab);

    if (contains(verbosity,VerbosityLevel::debug))
    	out << "CartesianGrid<" << dim << "> memory address = " << grid_ << endl;

    out << "Flat id = " << this->get_flat_index() << endl;
    out << "Tensor id = " << this->get_tensor_index() << endl;

    out << "Box intervals: " << endl;
    out.push(tab);
    for (int i = 0 ; i < dim ; ++i)
    {
        const auto coord_a = this->vertex(0)[i];
        const auto coord_b = this->vertex(pow(2,i))[i];
        out << "Direction["<< i << "] : [ " << coord_a << " , " << coord_b << " ]" << endl;
    }
    out.pop();



    out << "Vertices:" << endl;
    out.push(tab);
    for (Size i = 0 ; i < UnitElement<dim>::vertices_per_element ; ++i)
    {
        out << "Vertex[" << i << "] = " << this->vertex(i) << endl;
    }
    out.pop();

    out << "Center = " << this->center() << endl;

    std::string is_boundary = (this->is_boundary())?"true":"false";
    out << "On boundary = " << is_boundary << endl;

    if (this->is_boundary())
    {
        out << "Faces on boundary =";
        for (Size i = 0 ; i < UnitElement<dim>::faces_per_element ; ++i)
            if (this->is_boundary(i))
                out << " " << i ;

        out << endl;
    }

    out.pop();
}



IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid_element.inst>
