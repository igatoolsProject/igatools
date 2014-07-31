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

template <int dim_>
CartesianGridElement<dim_>::
CartesianGridElement(const std::shared_ptr<ContainerType> grid,
                     const Index index)
    :
    grid_(grid)
{
    this->reset_flat_tensor_indices(index);
}


template <int dim_>
auto
CartesianGridElement<dim_>::
get_grid() const -> const std::shared_ptr<ContainerType>
{
    return grid_;
}


template <int dim_>
inline
auto
CartesianGridElement<dim_>::
get_flat_index() const -> Index
{
    return flat_index_ ;
}


template <int dim_>
inline
auto
CartesianGridElement<dim_>::
get_tensor_index() const -> TensorIndex<dim>
{
    return tensor_index_ ;
}



template <int dim_>
void
CartesianGridElement<dim_>::
reset_flat_tensor_indices(const Index flat_index)
{
    Assert((flat_index == IteratorState::pass_the_end) ||
           ((flat_index >= 0) && (flat_index < grid_->get_num_all_elems())),
           ExcIndexRange(flat_index, 0, grid_->get_num_all_elems()));

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



template <int dim_>
void
CartesianGridElement<dim_>::
reset_flat_tensor_indices(const TensorIndex<dim> &tensor_index)
{
    tensor_index_= tensor_index;

    using Utils = MultiArrayUtils<dim>;
    flat_index_ = Utils::tensor_to_flat_index(
                      tensor_index_,
                      Utils::compute_weight(grid_->get_num_elements_dim()));

    Assert((flat_index_ == IteratorState::pass_the_end) ||
           ((flat_index_ >= 0) && (flat_index_ < grid_->get_num_active_elems())),
           ExcIndexRange(flat_index_, 0, grid_->get_num_active_elems()));
}



template <int dim_>
auto
CartesianGridElement<dim_>::
vertex(const int i) const -> Points<dim>
{
    Assert(i < UnitElement<dim>::vertices_per_element,
    ExcIndexRange(i,0,UnitElement<dim>::vertices_per_element));

    TensorIndex<dim> index = this->get_tensor_index();

    for (int j = 0; j < dim; ++j)
        index[j] += UnitElement<dim>::vertex_to_component[i][j];


    return this->get_grid()->get_knot_coordinates().cartesian_product(index);
};


template <int dim_>
auto
CartesianGridElement<dim_>::
get_coordinate_lengths() const -> array< Real, dim>
{
    const Points<dim> &p_origin = this->vertex(0);
    const Points<dim> &p_end = this->vertex(UnitElement<dim>::vertices_per_element-1);

    array<Real,dim> coord_length;
    for (int d = 0; d < dim ; ++d)
        coord_length[d] = p_end[d] - p_origin[d];

    return coord_length;
}


template <int dim_>
auto
CartesianGridElement<dim_>::
center() const -> Points< dim >
{
    Points<dim> center(vertex(0));
    center += vertex(UnitElement<dim>::opposite_vertex[0]);
    center *= 0.5;

    return (center) ;
}


template <int dim_>
Real
CartesianGridElement<dim_>::
get_measure(const TopologyId<dim> &topology_id) const
{
    const auto active_directions = topology_id.get_active_directions();

    Real result = 1.;
    for (const int &d : active_directions)
        result *= this->get_coordinate_lengths()[d];

    return result;
}



template <int dim_>
Real
CartesianGridElement<dim_>::
get_face_measure(const Index face_id) const
{
    return this->get_measure(FaceTopology<dim>(face_id));
}


template <int dim_>
bool
CartesianGridElement<dim_>::
is_point_inside(const Points< dim > &point) const
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

template <int dim_>
bool
CartesianGridElement<dim_>::
is_point_on_boundary(const Points< dim > &point) const
{
    int n_coords_inside = 0;
    int n_coords_on_boundary = 0;

    const auto &knots_directions = this->get_grid()->get_knot_coordinates();
    const auto &tensor_index = this->get_tensor_index();
    for (int j = 0; j < dim; ++j)
    {
        const auto &knots = knots_directions.get_data_direction(j) ;
        const Index id = tensor_index[j] ;

        const Real pt_coord = point(j) ;

        if (pt_coord > knots[id] && pt_coord < knots[id+1])
            n_coords_inside++;

        if (pt_coord == knots[id] || pt_coord == knots[id+1])
            n_coords_on_boundary++;
    }

    int n_coords_outside = dim - (n_coords_inside+n_coords_on_boundary);

    if (n_coords_on_boundary > 0 && n_coords_outside == 0)
        return true;
    else
        return false;
}


template <int dim_>
bool CartesianGridElement<dim_>::
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


template <int dim_>
bool
CartesianGridElement<dim_>::
is_boundary(const Index face_id) const
{
    const int const_direction = UnitElement<dim>::face_constant_direction[face_id];
    const int face_side = UnitElement<dim>::face_side[face_id];

    const auto element_id_dir = this->get_tensor_index()[const_direction] ;
    const auto num_elements_dir = this->get_grid()->get_num_elements_dim()(const_direction);

    return (element_id_dir == ((num_elements_dir-1) * face_side)) ;
}



template <int dim_>
auto
CartesianGridElement<dim_>::
transform_points_unit_to_reference(const vector<Points<dim>> &points_unit_domain) const ->
vector<Points<dim>>
{
    const int n_points = points_unit_domain.size();
    Assert(n_points > 0,ExcEmptyObject());


    const auto translate = this->vertex(0);
    const auto dilate    = this->get_coordinate_lengths();

    vector<Points<dim>> points_ref_domain(n_points);
    for (int ipt = 0 ; ipt < n_points ; ++ipt)
    {
        const auto &point_unit_domain = points_unit_domain[ipt];
        auto &point_ref_domain = points_ref_domain[ipt];
        for (int i = 0 ; i < dim ; ++i)
        {
            Assert(point_unit_domain[i] >= 0.0 && point_unit_domain[i] <= 1.0,
            ExcMessage("The coordinate of the point " + std::to_string(ipt) +
            " along the direction " + std::to_string(i) +
            " is not in the unit hypercube [0,1]^" + std::to_string(dim)));
            point_ref_domain[i] = point_unit_domain[i] * dilate[i] + translate[i];
        } // end loop i
    }

    return points_ref_domain;
}


template <int dim_>
auto
CartesianGridElement<dim_>::
transform_points_reference_to_unit(const vector<Points<dim>> &points_ref_domain) const ->
vector<Points<dim>>
{
    const int n_points = points_ref_domain.size();
    Assert(n_points > 0,ExcEmptyObject());


    const auto translate = this->vertex(0);
    const auto dilate    = this->get_coordinate_lengths();

    vector<Points<dim>> points_unit_domain(n_points);
//    LogStream out ;
//    using std::endl;
//    out << "CartesianGridElement::transform_points_reference_to_unit" << endl;
//    this->print_info(out);
//    out <<endl;
    for (int ipt = 0 ; ipt < n_points ; ++ipt)
    {
        const auto &point_ref_domain = points_ref_domain[ipt];
//        out << "point_ref_domain="<<point_ref_domain<<endl;
        Assert(this->is_point_inside(point_ref_domain) || this->is_point_on_boundary(point_ref_domain),
        ExcMessage("The point " + std::to_string(ipt) +
        " is outside this CartesianGridElement."));

        auto &point_unit_domain = points_unit_domain[ipt];
        for (int i = 0 ; i < dim ; ++i)
        {
//            out << "dilate["<<i<<"]=" <<dilate[i] << endl;
//            out << "translate["<<i<<"]=" <<translate[i] << endl;

            point_unit_domain[i] = (point_ref_domain[i] - translate[i]) / dilate[i] ;
        }
    }

    return points_unit_domain;
}

template <int dim_>
bool
CartesianGridElement<dim_>::
is_valid() const
{
    return (flat_index_ >= 0 && flat_index_ < grid_->get_num_active_elems())?true:false;
}


template <int dim_>
bool
CartesianGridElement<dim_>::
is_influence() const
{
    return grid_->influent_(flat_index_);
}

template <int dim_>
bool
CartesianGridElement<dim_>::
is_active() const
{
    return grid_->active_elems_(flat_index_);
}

template <int dim_>
void
CartesianGridElement<dim_>::
set_influence(const bool influence_flag)
{
    std::const_pointer_cast<CartesianGrid<dim>>(grid_)->
                                             influent_(flat_index_) = influence_flag;
}

template <int dim_>
void
CartesianGridElement<dim_>::
set_active(const bool active_flag)
{
    std::const_pointer_cast<CartesianGrid<dim>>(grid_)->
                                             active_elems_(flat_index_) = active_flag;
}

template <int dim_>
bool
CartesianGridElement<dim_>::
move(const TensorIndex<dim> &increment)
{
    tensor_index_ += increment;

    const auto n_elems = grid_->get_num_elements_dim();
    bool valid_tensor_index = true;
    for (int i = 0 ; i < dim ; ++i)
    {
        if (tensor_index_(i) < 0 || tensor_index_(i) >= n_elems(i))
        {
            valid_tensor_index = false;
            flat_index_ = IteratorState::invalid;
            break;
        }
    }

    if (valid_tensor_index)
    {
        using Utils = MultiArrayUtils<dim>;

        flat_index_ = Utils::tensor_to_flat_index(
                          tensor_index_,
                          Utils::compute_weight(n_elems));
    }

    return valid_tensor_index;
}

template <int dim_>
void
CartesianGridElement<dim_>::
print_info(LogStream &out, const VerbosityLevel verbosity) const
{
    using std::endl;

    const std::string tab = "   ";

    out << "CartesianGridElement<" << dim_ << "> info:" << endl;
    out.push(tab);

    if (contains(verbosity,VerbosityLevel::debug))
        out << "CartesianGrid<" << dim_ << "> memory address = " << grid_ << endl;

    out << "Flat id = " << this->get_flat_index() << endl;
    out << "Tensor id = " << this->get_tensor_index() << endl;

    out << "Box intervals: " << endl;
    out.push(tab);
    for (int i = 0 ; i < dim_ ; ++i)
    {
        const auto coord_a = this->vertex(0)[i];
        const auto coord_b = this->vertex(pow(2,i))[i];
        out << "Direction["<< i << "] : [ " << coord_a << " , " << coord_b << " ]" << endl;
    }
    out.pop();

    out << "Vertices:" << endl;
    out.push(tab);
    for (Size i = 0 ; i < UnitElement<dim_>::vertices_per_element ; ++i)
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
        for (Size i = 0 ; i < UnitElement<dim_>::faces_per_element ; ++i)
            if (this->is_boundary(i))
                out << " " << i ;

        out << endl;
    }

    out.pop();
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid_element.inst>
