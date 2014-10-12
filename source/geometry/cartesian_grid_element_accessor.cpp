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


#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/geometry/unit_element.h>
#include <algorithm>

using std::array;
using std::endl;

IGA_NAMESPACE_OPEN

template <int dim_>
CartesianGridElement<dim_>::
CartesianGridElement(const std::shared_ptr<ContainerType> grid,
                     const Index index)
    :
    grid_(grid)
{
    move_to(index);
}



template <int dim_>
CartesianGridElement<dim_>::
CartesianGridElement(const std::shared_ptr<ContainerType> grid,
                     const TensorIndex<dim> index)
    :
    grid_(grid)
{
    move_to(index);
}



template <int dim_>
CartesianGridElement<dim_>::
CartesianGridElement(const CartesianGridElement<dim_> &elem, const CopyPolicy &copy_policy)
{
    flat_index_   = elem.flat_index_;
    tensor_index_ = elem.tensor_index_;

    if (elem.local_cache_ != nullptr)
    {
        if (copy_policy == CopyPolicy::shallow)
        {
            local_cache_ = elem.local_cache_;
        }
        else
        {
            local_cache_ = std::shared_ptr<LocalCache>(new LocalCache(*elem.local_cache_));
        }
    }
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
move_to(const Index flat_index)
{
    Assert((flat_index == IteratorState::pass_the_end) ||
           ((flat_index >= 0) && (flat_index < grid_->get_num_all_elems())),
           ExcIndexRange(flat_index, 0, grid_->get_num_all_elems()));

    flat_index_ = flat_index ;

    if (flat_index_ != IteratorState::pass_the_end)
        tensor_index_ = grid_->flat_to_tensor(flat_index_);
    else
        tensor_index_.fill(IteratorState::pass_the_end);
}



template <int dim_>
void
CartesianGridElement<dim_>::
move_to(const TensorIndex<dim> &tensor_index)
{
    tensor_index_= tensor_index;
    flat_index_ = grid_->tensor_to_flat(tensor_index_);

    Assert((flat_index_ == IteratorState::pass_the_end) ||
           ((flat_index_ >= 0) && (flat_index_ < grid_->get_num_active_elems())),
           ExcIndexRange(flat_index_, 0, grid_->get_num_active_elems()));
}



template <int dim_>
bool
CartesianGridElement<dim_>::
jump(const TensorIndex<dim> &increment)
{
    tensor_index_ += increment;

    const auto n_elems = grid_->get_num_intervals();
    bool valid_tensor_index = true;
    for (int i = 0 ; i < dim ; ++i)
        if (tensor_index_[i] < 0 || tensor_index_[i] >= n_elems[i])
        {
            valid_tensor_index = false;
            flat_index_ = IteratorState::invalid;
            break;
        }

    if (valid_tensor_index)
        flat_index_ = grid_->tensor_to_flat(tensor_index_);

    return valid_tensor_index;
}


template <int dim_>
void
CartesianGridElement<dim_>::
operator++()
{
    const auto n_elem = this->grid_->get_num_all_elems();
    Index index = this->get_flat_index();
    do
    {
        ++index;
    }
    while (index<n_elem && (!this->grid_->active_elems_[index]));

    if (index >= n_elem)
        index = IteratorState::pass_the_end;

    this->move_to(index);
}



template <int dim_>
bool
CartesianGridElement<dim_>::
is_influence() const
{
    return grid_->marked_elems_[flat_index_];
}



template <int dim_>
bool
CartesianGridElement<dim_>::
is_active() const
{
    return grid_->active_elems_[flat_index_];
}



template <int dim_>
void
CartesianGridElement<dim_>::
set_influence(const bool influence_flag)
{
    std::const_pointer_cast<CartesianGrid<dim>>(grid_)->
                                             marked_elems_[flat_index_] = influence_flag;
}



template <int dim_>
void
CartesianGridElement<dim_>::
set_active(const bool active_flag)
{
    std::const_pointer_cast<CartesianGrid<dim> >(grid_)->
    active_elems_[flat_index_] = active_flag;
}


template <int dim_>
bool
CartesianGridElement<dim_>::
operator==(const CartesianGridElement<dim_> &elem) const
{
    Assert(this->get_grid() == elem.get_grid(), ExcMessage("Cannot compare elements on different grid."));
    return (this->get_flat_index() == elem.get_flat_index());
}



template <int dim_>
bool
CartesianGridElement<dim_>::
operator!=(const CartesianGridElement<dim_> &elem) const
{
    Assert(this->get_grid() == elem.get_grid(), ExcMessage("Cannot compare elements on different grid."));
    return (this->get_flat_index() != elem.get_flat_index());
}



template <int dim_>
CartesianGridElement<dim_> &
CartesianGridElement<dim_>::
operator=(const CartesianGridElement<dim_> &element)
{
    shallow_copy_from(element);
    return *this;
}


template <int dim_>
void
CartesianGridElement<dim_>::
copy_from(const CartesianGridElement<dim_> &elem,
          const CopyPolicy &copy_policy)
{
    // TODO (pauletti, Oct 9, 2014): assert(equal_grids);
    flat_index_   = elem.flat_index_;
    tensor_index_ = elem.tensor_index_;
    if (this != &elem)
    {
        if (copy_policy == CopyPolicy::deep)
        {
            Assert(elem.local_cache_ != nullptr, ExcNullPtr());
            local_cache_ = std::shared_ptr<LocalCache>(new LocalCache(*elem.local_cache_));
        }
        else if (copy_policy == CopyPolicy::shallow)
        {
            local_cache_ = elem.local_cache_;
        }
        else
        {
            Assert(false,ExcNotImplemented());
            AssertThrow(false,ExcNotImplemented());
        }
    }
}



template <int dim_>
void
CartesianGridElement<dim_>::
deep_copy_from(const CartesianGridElement<dim_> &elem)
{
    copy_from(elem,CopyPolicy::deep);
}



template <int dim_>
void
CartesianGridElement<dim_>::
shallow_copy_from(const CartesianGridElement<dim_> &elem)
{
    copy_from(elem,CopyPolicy::shallow);
}


//template <int dim_>
//void
//CartesianGridElement<dim_>::
//LengthCache::
//resize(const CartesianGrid<dim_> &grid)
//{
//    length_data_ = grid.get_element_lengths();
//
//    auto const size = length_data_.tensor_size();
//    length_.resize(size);
//    this->set_initialized(true);
//
//
//    for (int i = 0; i < dim_; ++i)
//        for (int j = 0; j < size(i); ++j)
//            length_.entry(i,j) = &length_data_.entry(i,j);
//
//    this->set_filled(true);
//}



template <int dim_>
auto
CartesianGridElement<dim_>::
vertex(const int i) const -> Point
{
    Assert(i < UnitElement<dim>::vertices_per_element,
           ExcIndexRange(i,0,UnitElement<dim>::vertices_per_element));

    TensorIndex<dim> index = this->get_tensor_index();

    for (int j = 0; j < dim; ++j)
        index[j] += UnitElement<dim>::vertex_to_component[i][j];


    return grid_->get_knot_coordinates().cartesian_product(index);
}



//template <int dim_>
//auto
//CartesianGridElement<dim_>::
//get_values_cache(const TopologyId<dim_> &topology_id) const -> const ValuesCache &
//{
//    Assert(topology_id.is_element() || topology_id.is_face(),
//           ExcMessage("Only element or face topology is allowed."));
//    Assert(local_cache_ != nullptr,ExcNullPtr());
//    if (topology_id.is_element())
//    {
//        return local_cache_->elem_values_;
//    }
//    else
//    {
//        // TODO (pauletti, Oct 9, 2014): reinsert this assertion
////        Assert(this->is_boundary(topology_id.get_id()),
////               ExcMessage("The requested face_id=" +
////                          std::to_string(topology_id.get_id()) +
////                          " is not a boundary for the element"));
//
//        return local_cache_->face_values_[topology_id.get_id()];
//    }
//}

//template <int dim_>
//auto
//CartesianGridElement<dim_>::
//get_values_cache(const TopologyId<dim_> &topology_id) -> ValuesCache &
//{
//    Assert(topology_id.is_element() || topology_id.is_face(),
//    ExcMessage("Only element or face topology is allowed."));
//    Assert(local_cache_ != nullptr,ExcNullPtr());
//    if (topology_id.is_element())
//    {
//        return local_cache_->elem_values_;
//    }
//    else
//    {
//        // TODO (pauletti, Oct 9, 2014): reinsert this assertion
////        Assert(this->is_boundary(topology_id.get_id()),
////        ExcMessage("The requested face_id=" +
////        std::to_string(topology_id.get_id()) +
////        " is not a boundary for the element"));
//
//        return local_cache_->face_values_[topology_id.get_id()];
//    }
//}


//template <int dim_>
//void
//CartesianGridElement<dim_>::
//init_cache(const ValueFlags flag,
//           const Quadrature<dim_> &quad)
//{
//    Assert((flag|admisible_flag) == admisible_flag,
//           ExcFillFlagNotSupported(admisible_flag, flag));
//    Assert(length_cache_.use_count() == 1, ExcCacheInUse(length_cache_.use_count()));
//
//
//    length_cache_->resize(*this->get_grid());
//
//    GridElemValueFlagsHandler elem_flags_handler(flag);
//    GridFaceValueFlagsHandler face_flags_handler(flag);
//
//    elem_values_.resize(elem_flags_handler, quad);
//
//    Index face_id = 0 ;
//    for (auto &face_value : face_values_)
//        face_value.resize(face_flags_handler, quad, face_id++);
//}



//template <int dim_>
//void
//CartesianGridElement<dim_>::
//init_cache(const ValueFlags flag)
//{
//    length_cache_->resize(*this->get_grid());
//
//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
//}

template <int dim_>
bool CartesianGridElement<dim_>::
is_boundary() const
{
    const auto num_elements_dim = this->get_grid()->get_num_intervals();

    const auto &element_index = this->get_tensor_index() ;

    for (int i = 0; i < dim; ++i)
    {
        if (element_index[i] == 0 or element_index[i] == num_elements_dim[i] - 1)
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
    const auto num_elements_dir = this->get_grid()->get_num_intervals()[const_direction];

    return (element_id_dir == ((num_elements_dir-1) * face_side)) ;
}



template <int dim_>
template <int k>
Real
CartesianGridElement<dim_>::
get_measure_(const int j) const
{
    const auto &cache = local_cache_->template get_value_cache<k>(j);
    Assert(cache.is_filled(), ExcMessage("Cache not filed."));
    Assert(cache.flags_handler_.measures_filled(), ExcMessage("Cache not filed."));

    return cache.measure_;
}



template <int dim_>
inline Real
CartesianGridElement<dim_>::
get_measure() const
{
    return get_measure_<0>(0);
}


template <int dim_>
inline Real
CartesianGridElement<dim_>::
get_face_measure(const Index face_id) const
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    return get_measure_<1>(face_id);
}



template <int dim_>
template <int k>
ValueVector<Real>
CartesianGridElement<dim_>::
get_w_measures_(const int j) const
{
    const auto &cache = local_cache_->template get_value_cache<k>(j);
    Assert(cache.is_filled(), ExcNotInitialized());
    Assert(cache.flags_handler_.measures_filled(), ExcNotInitialized());
    //Assert(cache.flags_handler_.weights_filled(), ExcNotInitialized());
    return (cache.measure_ * cache.unit_weights_);
}



template <int dim_>
ValueVector<Real>
CartesianGridElement<dim_>::
get_w_measures() const
{
    return get_w_measures_<0>(0);
}



template <int dim_>
ValueVector<Real>
CartesianGridElement<dim_>::
get_face_w_measures(const Index face_id) const
{
    return get_w_measures_<1>(face_id);
}



template <int dim_>
template <int k>
auto
CartesianGridElement<dim_>::
get_coordinate_lengths_(const int j) const -> const Point &
{
    const auto &cache = local_cache_->template get_value_cache<k>(j);
    Assert(cache.is_filled(), ExcNotInitialized());
    Assert(cache.flags_handler_.lengths_filled(), ExcNotInitialized());
    return cache.lengths_;
}


template <int dim_>
auto
CartesianGridElement<dim_>::
get_coordinate_lengths() const -> const Point &
{
    return get_coordinate_lengths_<0>(0);
}


template <int dim_>
template <int k>
auto
CartesianGridElement<dim_>::get_points_(const int j) const ->ValueVector<Point>
{
    const auto &cache =  local_cache_->template get_value_cache<k>(j);
    Assert(cache.flags_handler_.points_filled(), ExcNotInitialized());
    auto translate = vertex(0);
    auto dilate    = get_coordinate_lengths_<k>(j);

    auto ref_points = cache.unit_points_;
    ref_points.dilate_translate(dilate, translate);

    return ref_points.get_flat_cartesian_product();
}



template <int dim_>
auto
CartesianGridElement<dim_>::
get_points() const -> ValueVector<Point> const
{
    return get_points_<0>(0);
}



template <int dim_>
auto
CartesianGridElement<dim_>::
get_face_points(const Index face_id) const -> ValueVector<Point> const
{
    return get_points_<1>(face_id);
}


//template <int dim_>
//array< Real, dim_>
//CartesianGridElement<dim_>::
//get_coordinate_lengths() const
//{
//    Assert(length_cache_->is_filled(),ExcMessage("Cache not filled"));
//
//    const auto &tensor_index = this->get_tensor_index();
//
//    array<Real,dim_> coord_length;
//    for (int d = 0; d<dim_; d++)
//    {
//        const auto &length_d = length_cache_->length_.get_data_direction(d);
//        coord_length[d] = *(length_d[tensor_index[d]]);
//    }
//    return coord_length;
//}



template <int dim_>
void
CartesianGridElement<dim_>::
ValuesCache::
resize(const GridElemValueFlagsHandler &flags_handler,
       const Quadrature<dim> &quad)
{
    //const auto n_points_direction = quad.get_num_points_direction();
    //const Size n_points = n_points_direction.flat_size();

    flags_handler_ = flags_handler;

    if (flags_handler_.fill_points())
    {
        this->unit_points_ = quad.get_points();
        flags_handler_.set_points_filled(true);
    }

    if (flags_handler_.fill_w_measures())
    {
        this->unit_weights_ = quad.get_weights().get_flat_tensor_product();
    }
    else
    {
        this->unit_weights_.clear() ;
    }
    this->set_initialized(true);
}



template <int dim_>
void
CartesianGridElement<dim_>::
ValuesCache::print_info(LogStream &out) const
{
    out.begin_item("Fill flags:");
    flags_handler_.print_info(out);
    out.end_item();

    out << "Measure: " << measure_ << std::endl;
    out << "Lengths: " << lengths_ << std::endl;
    out.begin_item("Unit weights:");
    unit_weights_.print_info(out);
    out.end_item();

    out.begin_item("Unit points:");
    unit_points_.print_info(out);
    out.end_item();
}



template <int dim_>
void
CartesianGridElement<dim_>::
print_info(LogStream &out) const
{
    out << "Flat id = "   << flat_index_ << "    ";
    out << "Tensor id = " << tensor_index_ << endl;
}



template <int dim_>
void
CartesianGridElement<dim_>::
LocalCache::
print_info(LogStream &out) const
{
    out.begin_item("Element Cache:");
    std::get<0>(values_)[0].print_info(out);
    out.end_item();

    for (int i = 0 ; i < n_faces ; ++i)
    {
        out.begin_item("Face: "+ std::to_string(i) + " Cache:");
        std::get<1>(values_)[i].print_info(out);
        out.end_item();
    }
}



template <int dim_>
void
CartesianGridElement<dim_>::
print_cache_info(LogStream &out) const
{
    Assert(local_cache_ != nullptr, ExcNullPtr());
    local_cache_->print_info(out);
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid_element_accessor.inst>
