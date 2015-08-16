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
CartesianGridElement(const CartesianGridElement<dim> &elem, const CopyPolicy &copy_policy)
{
    grid_         = elem.grid_;
    tensor_index_ = elem.tensor_index_;

    if (elem.all_sub_elems_cache_ != nullptr)
    {
        if (copy_policy == CopyPolicy::shallow)
        {
            all_sub_elems_cache_ = elem.all_sub_elems_cache_;
        }
        else
        {
            all_sub_elems_cache_ = std::make_shared<CacheType>(*elem.all_sub_elems_cache_);
        }
    }
}


template <int dim>
std::shared_ptr<CartesianGridElement<dim> >
CartesianGridElement<dim>::
clone() const
{
//    auto elem = std::shared_ptr<CartesianGridElement<dim> >(
//                   new CartesianGridElement(*this,CopyPolicy::deep));
    auto elem = std::make_shared<CartesianGridElement<dim>>(*this,CopyPolicy::deep);
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
bool
CartesianGridElement<dim>::
is_property_true(const std::string &property) const
{
    const auto &elems_same_property = grid_->get_elements_id_same_property(property);
    return std::binary_search(elems_same_property.begin(),elems_same_property.end(),flat_index_);
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
            Assert(elem.all_sub_elems_cache_ != nullptr, ExcNullPtr());
            all_sub_elems_cache_ = std::make_shared<CacheType>(*elem.all_sub_elems_cache_);
        }
        else if (copy_policy == CopyPolicy::shallow)
        {
            all_sub_elems_cache_ = elem.all_sub_elems_cache_;
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
    for (const auto j : UnitElement<dim>::active_directions)
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

    auto &k_elem = UnitElement<dim>::template get_elem<k>(id);

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
    for (auto &id : UnitElement<dim>::template elems_ids<k>())
        if (is_boundary<k>(id))
            return true;

    return false;
}




template <int dim>
template <int k>
Real
CartesianGridElement<dim>::
get_measure(const int j) const
{
    const auto lengths = this->template get_side_lengths<k>(j);

    auto &k_elem = UnitElement<dim>::template get_elem<k>(j);

    Real measure = 1.0;
    for (const int active_dir : k_elem.active_directions)
        measure *= lengths[active_dir];

    return measure;
}





template <int dim>
template <int k>
ValueVector<Real>
CartesianGridElement<dim>::
get_w_measures(const int j) const
{
    return this->template get_values_from_cache<_W_Measure,k>(j);
}




template <int dim>
template <int k>
auto
CartesianGridElement<dim>::
get_side_lengths(const int j) const -> const Point
{
    Point lengths;
#if 0
    auto &k_elem = UnitElement<dim>::template get_elem<k>(j);

    for (const int const_dir :k_elem.constant_directions)
        lengths[const_dir] = 0.0;

    for (const int active_dir : k_elem.active_directions)
    {
        const auto &knots_active_dir = grid_->get_knot_coordinates(active_dir);
        const int j = tensor_index_[active_dir];
        lengths[active_dir] = knots_active_dir[j+1] - knots_active_dir[j];
    }
#endif
    for (const int active_dir : UnitElement<dim>::active_directions)
    {
        const auto &knots_active_dir = grid_->get_knot_coordinates(active_dir);
        const int j = tensor_index_[active_dir];
        lengths[active_dir] = knots_active_dir[j+1] - knots_active_dir[j];
    }

    return lengths;
}



template <int dim>
template <int k>
auto
CartesianGridElement<dim>::
get_points(const int j) const ->ValueVector<Point>
{
    return this->template get_values_from_cache<_Point,k>(j);
}




template <int dim>
auto
CartesianGridElement<dim>::
get_element_points() const -> ValueVector<Point>
{
    return this->template get_points<dim>(0);
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
CartesianGridElement<dim>::
print_cache_info(LogStream &out) const
{
//    Assert(all_sub_elems_cache_ != nullptr, ExcNullPtr());
    if (all_sub_elems_cache_)
        all_sub_elems_cache_->print_info(out);
    else
        out << "Cache not allocated." << std::endl;
}


template <int dim>
SafeSTLVector<std::string>
CartesianGridElement<dim>::
get_defined_properties() const
{
    SafeSTLVector<std::string> elem_properties;

    SafeSTLVector<std::string> grid_properties = grid_->properties_elements_id_.get_properties();
    for (const auto &property : grid_properties)
    {
        if (grid_->test_if_element_has_property(flat_index_, property))
            elem_properties.emplace_back(property);
    }
    return elem_properties;
}



template <int dim>
ValueFlags
CartesianGridElement<dim>::
get_valid_flags()
{
    return cacheutils::get_valid_flags_from_cache_type(CType());
}


#ifdef SERIALIZATION
template <int dim>
template<class Archive>
void
CartesianGridElement<dim>::
serialize(Archive &ar, const unsigned int version)
{
    auto non_const_grid = std::const_pointer_cast<CartesianGrid<dim>>(grid_);
    ar &boost::serialization::make_nvp("grid_",non_const_grid);
    grid_ = non_const_grid;
    Assert(grid_ != nullptr,ExcNullPtr());

    ar &boost::serialization::make_nvp("flat_index_",flat_index_);

    ar &boost::serialization::make_nvp("tensor_index_",tensor_index_);

    ar &boost::serialization::make_nvp("all_sub_elems_cache_",all_sub_elems_cache_);
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid_element.inst>
