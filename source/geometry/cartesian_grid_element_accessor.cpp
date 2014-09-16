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

IGA_NAMESPACE_OPEN

//TODO: inline the appropriate method and put in separate file




template <int dim_>
CartesianGridElementAccessor<dim_>::
CartesianGridElementAccessor(const std::shared_ptr<ContainerType> grid,
                             const Index index)
    :
    CartesianGridElement<dim>(grid, index)
{}



template <int dim_>
CartesianGridElementAccessor<dim_>::
CartesianGridElementAccessor(const std::shared_ptr<ContainerType> grid,
                             const TensorIndex<dim> index)
    :
    CartesianGridElement<dim>(grid, index)
{}

template <int dim_>
CartesianGridElementAccessor<dim_>::
CartesianGridElementAccessor(const CartesianGridElementAccessor<dim_> &elem, const CopyPolicy &copy_policy)
    :
    CartesianGridElement<dim_>(elem)
{
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
void
CartesianGridElementAccessor<dim_>::
copy_from(const CartesianGridElementAccessor<dim_> &elem,
          const CopyPolicy &copy_policy)
{
    parent_t::operator=(elem);
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
CartesianGridElementAccessor<dim_>::
deep_copy_from(const CartesianGridElementAccessor<dim_> &elem)
{
    this->copy_from(elem,CopyPolicy::deep);
}

template <int dim_>
void
CartesianGridElementAccessor<dim_>::
shallow_copy_from(const CartesianGridElementAccessor<dim_> &elem)
{
    this->copy_from(elem,CopyPolicy::shallow);
}


template <int dim_>
CartesianGridElementAccessor<dim_> &
CartesianGridElementAccessor<dim_>::
operator=(const CartesianGridElementAccessor<dim_> &element)
{
    this->shallow_copy_from(element);
    return (*this);
}


//template <int dim_>
//void
//CartesianGridElementAccessor<dim_>::
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
CartesianGridElementAccessor<dim_>::
get_values_cache(const TopologyId<dim_> &topology_id) const -> const ValuesCache &
{
    Assert(topology_id.is_element() || topology_id.is_face(),
           ExcMessage("Only element or face topology is allowed."));
    Assert(local_cache_ != nullptr,ExcNullPtr());
    if (topology_id.is_element())
    {
        return local_cache_->elem_values_;
    }
    else
    {
        Assert(this->is_boundary(topology_id.get_id()),
               ExcMessage("The requested face_id=" +
                          std::to_string(topology_id.get_id()) +
                          " is not a boundary for the element"));

        return local_cache_->face_values_[topology_id.get_id()];
    }
}

template <int dim_>
auto
CartesianGridElementAccessor<dim_>::
get_values_cache(const TopologyId<dim_> &topology_id) -> ValuesCache &
{
    Assert(topology_id.is_element() || topology_id.is_face(),
    ExcMessage("Only element or face topology is allowed."));
    Assert(local_cache_ != nullptr,ExcNullPtr());
    if (topology_id.is_element())
    {
        return local_cache_->elem_values_;
    }
    else
    {
        Assert(this->is_boundary(topology_id.get_id()),
        ExcMessage("The requested face_id=" +
        std::to_string(topology_id.get_id()) +
        " is not a boundary for the element"));

        return local_cache_->face_values_[topology_id.get_id()];
    }
}


//template <int dim_>
//void
//CartesianGridElementAccessor<dim_>::
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
//CartesianGridElementAccessor<dim_>::
//init_cache(const ValueFlags flag)
//{
//    length_cache_->resize(*this->get_grid());
//
//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
//}


template <int dim_>
inline Real
CartesianGridElementAccessor<dim_>::
get_measure(const TopologyId<dim_> &topology_id) const
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled(), ExcMessage("Cache not filed."));
    Assert(cache.flags_handler_.measures_filled(), ExcMessage("Cache not filed."));

    return cache.measure_;
}


template <int dim_>
inline Real
CartesianGridElementAccessor<dim_>::
get_face_measure(const Index face_id) const
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    return this->get_measure(FaceTopology<dim_>(face_id));
}





//template <int dim_>
//void
//CartesianGridElementAccessor<dim_>::
//init_face_cache(const Index face_id,
//                const ValueFlags flag,
//                const Quadrature<dim_-1> &quad)
//{
//    Assert(false, ExcNotImplemented());
//}



//template <int dim_>
//void
//CartesianGridElementAccessor<dim_>::
//fill_cache(const TopologyId<dim_> &topology_id)
//{
//    auto &cache = get_values_cache(topology_id);
//
//    cache.fill(CartesianGridElement<dim_>::get_measure(topology_id));
//
//    cache.set_filled(true);
//}



//template <int dim_>
//void
//CartesianGridElementAccessor<dim_>::
//fill_face_cache(const Index face_id)
//{
//    fill_cache(FaceTopology<dim_>(face_id));
//}




template <int dim_>
ValueVector<Real> const &
CartesianGridElementAccessor<dim_>::
get_w_measures(const TopologyId<dim_> &topology_id) const
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled(), ExcNotInitialized());
    Assert(cache.flags_handler_.w_measures_filled(), ExcNotInitialized());
    return cache.w_measure_;
}

template <int dim_>
ValueVector<Real> const &
CartesianGridElementAccessor<dim_>::
get_face_w_measures(const Index face_id) const
{
    return this->get_w_measures(FaceTopology<dim_>(face_id));
}


template <int dim_>
auto
CartesianGridElementAccessor<dim_>::
get_points(const TopologyId<dim_> &topology_id) const -> ValueVector<Points<dim>> const
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.flags_handler_.points_filled(), ExcNotInitialized());
    auto translate = this->vertex(0);
    auto dilate    = this->get_coordinate_lengths();

    auto ref_points = cache.unit_points_;
    ref_points.dilate_translate(dilate, translate);

    return ref_points.get_flat_cartesian_product();
}


template <int dim_>
auto
CartesianGridElementAccessor<dim_>::
get_face_points(const Index face_id) const -> ValueVector<Points<dim>> const
{
    return this->get_points(FaceTopology<dim_>(face_id));
}


//template <int dim_>
//array< Real, dim_>
//CartesianGridElementAccessor<dim_>::
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
CartesianGridElementAccessor<dim_>::
ValuesCache::
resize(const GridElemValueFlagsHandler &flags_handler,
       const Quadrature<dim_> &quad)
{
    const auto n_points_direction = quad.get_num_points_direction();
    const Size n_points = n_points_direction.flat_size();

    flags_handler_ = flags_handler;

    if (flags_handler_.fill_points())
    {
        this->unit_points_ = quad.get_points();
        flags_handler_.set_points_filled(true);
    }

    if (flags_handler_.fill_w_measures())
    {
        if (this->w_measure_.size() != n_points)
            this->w_measure_.resize(n_points);

        this->unit_weights_ = quad.get_weights().get_flat_tensor_product();
    }
    else
    {
        this->w_measure_.clear() ;
        this->unit_weights_.clear() ;
    }
    this->set_initialized(true);
}


template <int dim_>
void
CartesianGridElementAccessor<dim_>::
ValuesCache::
fill(const Real measure)
{
    if (flags_handler_.fill_measures())
    {
        this->measure_ = measure;

        flags_handler_.set_measures_filled(true);
    }

    if (flags_handler_.fill_w_measures())
    {
        Assert(flags_handler_.measures_filled(),ExcCacheNotFilled());

        w_measure_ = measure_ * unit_weights_;

        flags_handler_.set_w_measures_filled(true);
    }
}


template <int dim_>
void
CartesianGridElementAccessor<dim_>::
FaceValuesCache::
resize(const GridFaceValueFlagsHandler &flags_handler,const Quadrature<dim_> &quad, const Index face_id)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    ValuesCache::resize(flags_handler,quad.collapse_to_face(face_id));
}



template <int dim_>
void
CartesianGridElementAccessor<dim_>::
FaceValuesCache::
resize(const GridFaceValueFlagsHandler &flags_handler,const Quadrature<dim_-1> &quad1, const Index face_id)
{
    Assert(false, ExcNotImplemented());
}


template <int dim_>
void
CartesianGridElementAccessor<dim_>::
print_info(LogStream &out, const VerbosityLevel verbosity) const
{
    CartesianGridElement<dim_>::print_info(out,verbosity);

    if (local_cache_ != nullptr)
        local_cache_->print_info(out);
}



template <int dim_>
void
CartesianGridElementAccessor<dim_>::
LocalCache::
print_info(LogStream &out) const
{
    out.begin_item("Element Cache:");
    elem_values_.print_info(out);
    out.end_item();

    for (int i = 0 ; i < n_faces ; ++i)
    {
        out.begin_item("Face: "+ std::to_string(i) + " Cache:");
        face_values_[i].print_info(out);
        out.end_item();
    }
}

template <int dim_>
void
CartesianGridElementAccessor<dim_>::
print_cache_info(LogStream &out) const
{
    Assert(local_cache_ != nullptr, ExcNullPtr());
    local_cache_->print_info(out);
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid_element_accessor.inst>
