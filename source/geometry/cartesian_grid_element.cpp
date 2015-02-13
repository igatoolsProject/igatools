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


#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/geometry/unit_element.h>
#include <algorithm>

using std::array;
using std::endl;

IGA_NAMESPACE_OPEN

template <int dim>
CartesianGridElement<dim>::
CartesianGridElement(const std::shared_ptr<ContainerType> grid,
                     const Index index)
    :
    grid_(grid)
{
    Assert(grid_ != nullptr,ExcNullPtr());
    move_to(index);
}



template <int dim>
CartesianGridElement<dim>::
CartesianGridElement(const std::shared_ptr<ContainerType> grid,
                     const TensorIndex<dim> index)
    :
    CartesianGridElement(grid,grid->tensor_to_flat(index))
{}



template <int dim>
CartesianGridElement<dim>::
CartesianGridElement(const CartesianGridElement<dim> &elem, const CopyPolicy &copy_policy)
{
    grid_         = elem.grid_;
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


template <int dim>
std::shared_ptr<CartesianGridElement<dim> >
CartesianGridElement<dim>::
clone() const
{
    auto elem = std::shared_ptr<CartesianGridElement<dim> >(
                    new CartesianGridElement(*this,CopyPolicy::deep));
    Assert(elem != nullptr, ExcNullPtr());
    return elem;
}

template <int dim>
auto
CartesianGridElement<dim>::
get_grid() const -> const std::shared_ptr<const CartesianGrid<dim> >
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



template <int dim>
void
CartesianGridElement<dim>::
move_to(const TensorIndex<dim> &tensor_index)
{
    move_to(grid_->tensor_to_flat(tensor_index));
}



template <int dim>
bool
CartesianGridElement<dim>::
jump(const TensorIndex<dim> &increment)
{
    tensor_index_ += increment;

    const auto n_elems = grid_->get_num_intervals();
    bool valid_tensor_index = true;
    for (const auto i : Topology::active_directions)
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

#if 0
template <int dim>
void
CartesianGridElement<dim>::
operator++()
{
    const auto &active_elems = grid_->get_elements_id_same_property(
                                   CartesianGrid<dim>::ElementProperty::active);
    const auto elem_begin = active_elems.begin();
    const auto elem_end  = active_elems.end();

    Index index = this->get_flat_index();
    auto elem = std::find(elem_begin,elem_end,index);
    auto elem_next = ++elem;
    if (elem_next == elem_end)
        index = IteratorState::pass_the_end;
    else
        index = *elem_next;
    this->move_to(index);
}
#endif

template <int dim>
bool
CartesianGridElement<dim>::
is_property_true(const ElementProperty &property) const
{
    const auto &elems_same_property = grid_->get_elements_id_same_property(property);
    return std::binary_search(elems_same_property.begin(),elems_same_property.end(),flat_index_);
}



template <int dim>
bool
CartesianGridElement<dim>::
is_influence() const
{
    return is_property_true(ElementProperty::influence);
}



template <int dim>
bool
CartesianGridElement<dim>::
is_active() const
{
    return is_property_true(ElementProperty::active);
}



template <int dim>
void
CartesianGridElement<dim>::
set_influence(const bool influence_flag)
{
    using Grid = CartesianGrid<dim>;
    std::const_pointer_cast<Grid>(grid_)->set_element_property(
        ElementProperty::influence,
        flat_index_,
        influence_flag);
}





template <int dim>
bool
CartesianGridElement<dim>::
operator ==(const CartesianGridElement<dim> &elem) const
{
    Assert(this->get_grid() == elem.get_grid(), ExcMessage("Cannot compare elements on different grid."));
    return (this->get_flat_index() == elem.get_flat_index());
}



template <int dim>
bool
CartesianGridElement<dim>::
operator !=(const CartesianGridElement<dim> &elem) const
{
    Assert(this->get_grid() == elem.get_grid(), ExcMessage("Cannot compare elements on different grid."));
    return (this->get_flat_index() != elem.get_flat_index());
}

template <int dim>
bool
CartesianGridElement<dim>::
operator <(const CartesianGridElement<dim> &elem) const
{
    Assert(this->get_grid() == elem.get_grid(), ExcMessage("Cannot compare elements on different grid."));
    return (this->get_flat_index() < elem.get_flat_index());
}

template <int dim>
bool
CartesianGridElement<dim>::
operator >(const CartesianGridElement<dim> &elem) const
{
    Assert(this->get_grid() == elem.get_grid(), ExcMessage("Cannot compare elements on different grid."));
    return (this->get_flat_index() > elem.get_flat_index());
}


template <int dim>
CartesianGridElement<dim> &
CartesianGridElement<dim>::
operator=(const CartesianGridElement<dim> &element)
{
    shallow_copy_from(element);
    return *this;
}


template <int dim>
void
CartesianGridElement<dim>::
copy_from(const CartesianGridElement<dim> &elem,
          const CopyPolicy &copy_policy)
{
    Assert(this->get_grid() == elem.get_grid(), ExcMessage("Cannot copy from an element on different grid."));

    if (this != &elem)
    {
        flat_index_   = elem.flat_index_;
        tensor_index_ = elem.tensor_index_;

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



template <int dim>
void
CartesianGridElement<dim>::
deep_copy_from(const CartesianGridElement<dim> &elem)
{
    copy_from(elem,CopyPolicy::deep);
}



template <int dim>
void
CartesianGridElement<dim>::
shallow_copy_from(const CartesianGridElement<dim> &elem)
{
    copy_from(elem,CopyPolicy::shallow);
}



template <int dim>
auto
CartesianGridElement<dim>::
vertex(const int i) const -> Point
{
    Assert(i < UnitElement<dim>::sub_elements_size[0],
           ExcIndexRange(i,0, UnitElement<dim>::sub_elements_size[0]));

    TensorIndex<dim> index = this->get_tensor_index();

    auto all_elems = UnitElement<dim>::all_elems;
    for (const auto j : Topology::active_directions)
    {
        auto vertex = std::get<0>(all_elems)[i];
        index[j] += vertex.constant_values[j];
    }

    return grid_->get_knot_coordinates().cartesian_product(index);
}



template <int dim>
template <int k>
bool CartesianGridElement<dim>::
is_boundary(const Index id) const
{
    const auto &n_elem = this->get_grid()->get_num_intervals();
    const auto &index = this->get_tensor_index();

    auto &k_elem = Topology::template get_elem<k>(id);

    for (int i = 0; i < dim-k; ++i)
    {
        auto dir = k_elem.constant_directions[i];
        auto val = k_elem.constant_values[i];
        if (((index[dir] == 0)               && (val == 0)) ||
            ((index[dir] == n_elem[dir] - 1) && (val == 1)))
            return true;
    }

    return false;
}



template <int dim>
template <int k>
bool
CartesianGridElement<dim>::
is_boundary() const
{
    for (auto &id : Topology::template elems_ids<k>())
    {
        auto res = is_boundary<k>(id);
        if (res)
            return res;
    }
    return false;
}




template <int dim>
template <int k>
Real
CartesianGridElement<dim>::
get_measure(const int j) const
{
    const auto &cache = local_cache_->template get_value_cache<k>(j);
    Assert(cache.is_filled(), ExcMessage("Cache not filed."));
    Assert(cache.flags_handler_.measures_filled(), ExcMessage("Cache not filed."));

    return cache.measure_;
}





template <int dim>
template <int k>
ValueVector<Real>
CartesianGridElement<dim>::
get_w_measures(const int j) const
{
    const auto &cache = local_cache_->template get_value_cache<k>(j);
    Assert(cache.is_filled(), ExcNotInitialized());
    Assert(cache.flags_handler_.measures_filled(), ExcNotInitialized());
    //Assert(cache.flags_handler_.weights_filled(), ExcNotInitialized());
    return (cache.measure_ * cache.unit_weights_);
}




template <int dim>
template <int k>
auto
CartesianGridElement<dim>::
get_coordinate_lengths(const int j) const -> const Point &
{
    const auto &cache = local_cache_->template get_value_cache<k>(j);
    Assert(cache.is_filled(), ExcNotInitialized());
    Assert(cache.flags_handler_.lengths_filled(), ExcNotInitialized());
    return cache.lengths_;
}



template <int dim>
template <int k>
auto
CartesianGridElement<dim>::
get_points(const int j) const ->ValueVector<Point>
{
    const auto &cache = local_cache_->template get_value_cache<k>(j);
    Assert(cache.flags_handler_.points_filled(), ExcNotInitialized());
    auto translate = vertex(0);
    auto dilate    = get_coordinate_lengths<k>(j);


    const int n_pts = cache.unit_points_.get_num_points();
    ValueVector<Point> ref_points(n_pts);

    for (int ipt = 0 ; ipt < n_pts ; ++ipt)
    {
        const auto &unit_pt = cache.unit_points_[ipt];
        auto &ref_pt = ref_points[ipt];

        for (const auto dir : Topology::active_directions)
            ref_pt[dir] = unit_pt[dir] * dilate[dir] + translate[dir];
    }
    return ref_points;
}





template <int dim>
void
CartesianGridElement<dim>::
ValuesCache::
resize(const GridFlags &flags_handler,
       const EvaluationPoints<dim> &quad)
{
    flags_handler_ = flags_handler;

    if (flags_handler_.fill_points())
    {
        this->unit_points_ = quad.get_points();
        flags_handler_.set_points_filled(true);
    }

    if (flags_handler_.fill_w_measures())
    {
        this->unit_weights_ = quad.get_weights();
    }
    else
    {
        this->unit_weights_.clear() ;
    }
    this->set_initialized(true);
}



template <int dim>
void
CartesianGridElement<dim>::
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



template <int dim>
void
CartesianGridElement<dim>::
print_info(LogStream &out) const
{
    out << "Flat id = "   << flat_index_ << "    ";
    out << "Tensor id = " << tensor_index_ << endl;
}



template <int dim>
void
CartesianGridElement<dim>::LocalCache::
print_info(LogStream &out) const
{
    cacheutils::print_caches(values_, out);
//    out.begin_item("Element Cache:");
//    std::get<dim>(values_)[0].print_info(out);
//    out.end_item();

//    for (int i = 0 ; i < UnitElement<dim>::template num_elem<dim==0? 0 : dim-1>() ; ++i)
//    {
//        out.begin_item("Face: "+ std::to_string(i) + " Cache:");
//        std::get<1>(values_)[i].print_info(out);
//        out.end_item();
//    }
}



template <int dim>
void
CartesianGridElement<dim>::
print_cache_info(LogStream &out) const
{
    Assert(local_cache_ != nullptr, ExcNullPtr());
    local_cache_->print_info(out);
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid_element.inst>
