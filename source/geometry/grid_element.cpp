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


#include <igatools/geometry/grid_element.h>
#include <igatools/geometry/unit_element.h>
#include <algorithm>

IGA_NAMESPACE_OPEN

template <int dim, class ContainerType_>
GridElementBase<dim, ContainerType_>::
GridElementBase(const std::shared_ptr<ContainerType> grid,
                const ListIt &index,
                const PropId &prop)
    :
    grid_(grid),
    property_(prop),
    index_it_(index)
{}



template <int dim, class ContainerType_>
GridElementBase<dim, ContainerType_>::
GridElementBase(const self_t &elem,
                const CopyPolicy &copy_policy)
    :
    grid_(elem.grid_),
    property_(elem.property_),
    index_it_(elem.index_it_)
{
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



//template <int dim, class ContainerType_>
//std::shared_ptr<self_t >
//GridElementBase<dim, ContainerType_>::
//clone() const
//{
//    auto elem = std::make_shared<self_t>(*this,CopyPolicy::deep);
//    Assert(elem != nullptr, ExcNullPtr());
//    return elem;
//}



template <int dim, class ContainerType_>
auto
GridElementBase<dim, ContainerType_>::
get_grid() const -> const std::shared_ptr<const ContainerType>
{
    return grid_;
}



template <int dim, class ContainerType_>
auto
GridElementBase<dim, ContainerType_>::
get_index() const ->  const IndexType &
{
    return *index_it_;
}



//template <int dim, class ContainerType_>
//void
//GridElementBase<dim, ContainerType_>::
//move_to(const Index  &flat_index)
//{
//    Assert((flat_index == IteratorState::pass_the_end) ||
//           ((flat_index >= 0) && (flat_index < grid_->get_num_all_elems())),
//           ExcIndexRange(flat_index, 0, grid_->get_num_all_elems()));
//
//    flat_index_ = flat_index ;
//
//    if (flat_index_ != IteratorState::pass_the_end)
//        tensor_index_ = grid_->flat_to_tensor(flat_index_);
//    else
//        tensor_index_.fill(IteratorState::pass_the_end);
//}



template <int dim, class ContainerType_>
bool
GridElementBase<dim, ContainerType_>::
has_property(const PropId &prop) const
{
    const auto &list = grid_->elem_properties_[prop];
    return std::binary_search(list.begin(), list.end(), get_index());
}



template <int dim, class ContainerType_>
bool
GridElementBase<dim, ContainerType_>::
operator ==(const self_t &elem) const
{
    Assert(this->get_grid() == elem.get_grid(),
           ExcMessage("Cannot compare elements on different grid."));
    return (this->get_index() == elem.get_index());
}



template <int dim, class ContainerType_>
bool
GridElementBase<dim, ContainerType_>::
operator !=(const self_t &elem) const
{
    Assert(this->get_grid() == elem.get_grid(),
           ExcMessage("Cannot compare elements on different grid."));
    return (this->get_index() != elem.get_index());
}

template <int dim, class ContainerType_>
bool
GridElementBase<dim, ContainerType_>::
operator <(const self_t &elem) const
{
    Assert(this->get_grid() == elem.get_grid(),
           ExcMessage("Cannot compare elements on different grid."));
    return (this->get_index() < elem.get_index());
}

template <int dim, class ContainerType_>
bool
GridElementBase<dim, ContainerType_>::
operator >(const self_t &elem) const
{
    Assert(this->get_grid() == elem.get_grid(),
           ExcMessage("Cannot compare elements on different grid."));
    return (this->get_index() > elem.get_index());
}


template <int dim, class ContainerType_>
auto
GridElementBase<dim, ContainerType_>::
operator=(const self_t &element) -> self_t &
{
    shallow_copy_from(element);
    return *this;
}


template <int dim, class ContainerType_>
void
GridElementBase<dim, ContainerType_>::
copy_from(const self_t &elem, const CopyPolicy &copy_policy)
{
    Assert(this->get_grid() == elem.get_grid(),
           ExcMessage("Cannot copy from an element on different grid."));

    if (this != &elem)
    {
        index_it_   = elem.index_it_;
        property_ = elem.property_;

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



template <int dim, class ContainerType_>
void
GridElementBase<dim, ContainerType_>::
deep_copy_from(const self_t &elem)
{
    copy_from(elem,CopyPolicy::deep);
}



template <int dim, class ContainerType_>
void
GridElementBase<dim, ContainerType_>::
shallow_copy_from(const self_t &elem)
{
    copy_from(elem,CopyPolicy::shallow);
}



template <int dim, class ContainerType_>
auto
GridElementBase<dim, ContainerType_>::
vertex(const int i) const -> Point
{
    Assert(i < UnitElement<dim>::sub_elements_size[0],
           ExcIndexRange(i,0, UnitElement<dim>::sub_elements_size[0]));

    TensorIndex<dim> index = this->get_index();

    auto all_elems = UnitElement<dim>::all_elems;
    const auto &vertex = std::get<0>(all_elems)[i];

    for (const auto j : UnitElement<dim>::active_directions)
    {
        index[j] += vertex.constant_values[j];
    }

    return grid_->knot_coordinates_.cartesian_product(index);
}



template <int dim, class ContainerType_>
template <int k>
bool GridElementBase<dim, ContainerType_>::
is_boundary(const Index id) const
{
    const auto &n_elem = this->get_grid()->get_num_intervals();
    const auto &index = this->get_index();

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



template <int dim, class ContainerType_>
template <int k>
bool
GridElementBase<dim, ContainerType_>::
is_boundary() const
{
    for (auto &id : UnitElement<dim>::template elems_ids<k>())
        if (is_boundary<k>(id))
            return true;

    return false;
}




template <int dim, class ContainerType_>
template <int sdim>
Real
GridElementBase<dim, ContainerType_>::
get_measure(const int s_id) const
{
    const auto lengths = get_side_lengths<sdim>(s_id);

    //  auto &k_elem = UnitElement<dim>::template get_elem<k>(j);

    Real measure = 1.0;
    for (int i=0; i<sdim; ++i)
        measure *= lengths[i];

    return measure;
}




template <int dim, class ContainerType_>
template <int k>
ValueVector<Real>
GridElementBase<dim, ContainerType_>::
get_w_measures(const int j) const
{
    return this->template get_values_from_cache<_W_Measure,k>(j);
}



template <int dim, class ContainerType_>
template <int sdim>
auto
GridElementBase<dim, ContainerType_>::
get_side_lengths(const int sid) const -> const Points<sdim>
{
    Points<sdim> lengths;

    auto &s_elem = UnitElement<dim>::template get_elem<sdim>(sid);

    int i=0;
    for (const int active_dir : s_elem.active_directions)
    {
        const auto &knots_active_dir = grid_->get_knot_coordinates(active_dir);
        const int j = get_index()[active_dir];
        lengths[i] = knots_active_dir[j+1] - knots_active_dir[j];
        ++i;
    }

    return lengths;
}



template <int dim, class ContainerType_>
template <int k>
auto
GridElementBase<dim, ContainerType_>::
get_points(const int j) const ->ValueVector<Point>
{
    return this->template get_values_from_cache<_Point,k>(j);
}




template <int dim, class ContainerType_>
const std::string  GridElementBase<dim, ContainerType_>::_Point::name = "Element Quadrature Points";

template <int dim, class ContainerType_>
const std::string GridElementBase<dim, ContainerType_>::_W_Measure::name = "Element Quadrature Weights";



template <int dim, class ContainerType_>
auto
GridElementBase<dim, ContainerType_>::
get_element_points() const -> ValueVector<Point>
{
    return this->template get_points<dim>(0);
}


template <int dim, class ContainerType_>
void
GridElementBase<dim, ContainerType_>::
print_info(LogStream &out) const
{
    out.begin_item("Property: ");
    out << property_ << std::endl;
    out.end_item();
    out.begin_item("Index:");
    index_it_->print_info(out);
    out.end_item();
}



template <int dim, class ContainerType_>
void
GridElementBase<dim, ContainerType_>::
print_cache_info(LogStream &out) const
{
    if (all_sub_elems_cache_)
        all_sub_elems_cache_->print_info(out);
    else
        out << "Cache not allocated." << std::endl;
}



#if 0
template <int dim, class ContainerType_>
SafeSTLVector<std::string>
GridElementBase<dim, ContainerType_>::
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

#endif

#ifdef SERIALIZATION
template <int dim, class ContainerType_>
template<class Archive>
void
GridElementBase<dim, ContainerType_>::
serialize(Archive &ar, const unsigned int version)
{
    using namespace boost::serialization;
    auto non_const_grid = std::const_pointer_cast<CartesianGrid<dim>>(grid_);
    ar &make_nvp("grid_",non_const_grid);
    grid_ = non_const_grid;
    Assert(grid_ != nullptr, ExcNullPtr());

    ar &make_nvp("property_", property_);
    Assert(false, ExcNotImplemented());
    //ar &make_nvp("index_it_", index_it_);
    ar &make_nvp("quad_list_", quad_list_);
    ar &make_nvp("all_sub_elems_cache_",all_sub_elems_cache_);
}
#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_element.inst>
