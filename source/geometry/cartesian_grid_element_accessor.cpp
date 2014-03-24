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
using std::vector;

IGA_NAMESPACE_OPEN

//TODO: inline the appropriate method and put in separate file




template <int dim_>
CartesianGridElementAccessor<dim_>::
CartesianGridElementAccessor(
    const CartesianGrid<dim_> &patch,
    const Index index)
    :
    CartesianGridElement<dim>(patch,index),
    length_cache_ {new LengthCache}
{}



template <int dim_>
bool
CartesianGridElementAccessor<dim_>::
operator== (const CartesianGridElementAccessor<dim_> &a) const
{
    Assert(this->get_grid() == a.get_grid(), ExcMessage("Cannot Compare Iterators."));
    return (this->get_flat_index() == a.get_flat_index());
}



template <int dim_>
bool
CartesianGridElementAccessor<dim_>::
operator!= (const CartesianGridElementAccessor<dim_> &a) const
{
    Assert(this->get_grid() == a.get_grid(), ExcMessage("Cannot Compare Iterators."));
    return (this->get_flat_index() != a.get_flat_index());
}




template <int dim_>
void
CartesianGridElementAccessor<dim_>::
operator ++ ()
{
    Index index = this->get_flat_index();
    ++index;
    if (index >= this->get_grid()->get_num_elements())
        index = IteratorState::pass_the_end;

    this->reset_flat_tensor_indices(index);
}




template <int dim_>
void
CartesianGridElementAccessor<dim_>::
LengthCache::
reset(const CartesianGrid<dim_> &grid)
{
    length_data_ = grid.get_element_lengths();

    auto const size = length_data_.tensor_size();
    length_.resize(size);
    this->set_initialized(true);


    for (int i = 0; i < dim_; ++i)
        for (int j = 0; j < size(i); ++j)
            length_.entry(i,j) = &length_data_.entry(i,j);

    this->set_filled(true);
}



template <int dim_>
auto
CartesianGridElementAccessor<dim_>::
get_values_cache(const TopologyId<dim_> &topology_id) const -> const ValuesCache &
{
    Assert(topology_id.is_element() || topology_id.is_face(),
           ExcMessage("Only element or face topology is allowed."));
    if (topology_id.is_element())
    {
        return elem_values_;
    }
    else
    {
        return face_values_[topology_id.get_id()];
    }
}



template <int dim_>
void
CartesianGridElementAccessor<dim_>::
init_values(const ValueFlags flag,
            const Quadrature<dim_> &quad)
{
    Assert((flag|admisible_flag) == admisible_flag,
           ExcFillFlagNotSupported(admisible_flag, flag));

    length_cache_->reset(*this->get_grid());


    GridElemValueFlagsHandler elem_flags_handler(flag);
    GridFaceValueFlagsHandler face_flags_handler(flag);

    elem_values_.reset(elem_flags_handler,quad);

    Index face_id = 0 ;
    for (auto& face_value : face_values_)
        face_value.reset(face_flags_handler, quad, face_id++);
}

template <int dim_>
void
CartesianGridElementAccessor<dim_>::
init_values(const ValueFlags flag)
{
    length_cache_->reset(*this->get_grid());

    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}


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





template <int dim_>
void
CartesianGridElementAccessor<dim_>::
init_face_values(const Index face_id,
                 const ValueFlags flag,
                 const Quadrature<dim_-1> &quad)
{
    Assert(false, ExcNotImplemented());
}



template <int dim_>
void
CartesianGridElementAccessor<dim_>::
fill_values()
{
    elem_values_.fill(CartesianGridElement<dim_>::get_measure());

    elem_values_.set_filled(true);
}



template <int dim_>
void
CartesianGridElementAccessor<dim_>::
fill_face_values(const Index face_id)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    auto &face_value = face_values_[face_id] ;

    face_value.fill(CartesianGridElement<dim_>::get_face_measure(face_id));

    face_value.set_filled(true);
}




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
get_points(const TopologyId<dim_> &topology_id) const -> vector<Point<dim>> const
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
get_face_points(const Index face_id) const -> vector<Point<dim>> const
{
	return this->get_points(FaceTopology<dim_>(face_id));
}


template <int dim_>
array< Real, dim_>
CartesianGridElementAccessor<dim_>::
get_coordinate_lengths() const
{
    Assert(length_cache_->is_filled(),ExcMessage("Cache not filled"));

    const auto &tensor_index = this->get_tensor_index();

    array<Real,dim_> coord_length;
    for (int d = 0; d<dim_; d++)
    {
        const auto &length_d = length_cache_->length_.get_data_direction(d);
        coord_length[d] = *(length_d[tensor_index[d]]);
    }
    return coord_length;
}


template <int dim_>
void
CartesianGridElementAccessor<dim_>::
ValuesCache::
reset(const GridElemValueFlagsHandler &flags_handler,const Quadrature<dim_> &quad)
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
ElementValuesCache::
reset(const GridElemValueFlagsHandler &flags_handler,const Quadrature<dim_> &quad)
{
    ValuesCache::reset(flags_handler,quad);
    this->set_initialized(true);
}



template <int dim_>
void
CartesianGridElementAccessor<dim_>::
FaceValuesCache::
reset(const GridFaceValueFlagsHandler &flags_handler,const Quadrature<dim_> &quad1, const Index face_id)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    const auto quad = quad1.collapse_to_face(face_id);
    ValuesCache::reset(flags_handler,quad);
    this->set_initialized(true);
}



template <int dim_>
void
CartesianGridElementAccessor<dim_>::
FaceValuesCache::
reset(const GridFaceValueFlagsHandler &flags_handler,const Quadrature<dim_-1> &quad1, const Index face_id)
{
    Assert(false, ExcNotImplemented());
}


template <int dim_>
void
CartesianGridElementAccessor<dim_>::
print_info(LogStream &out, const VerbosityLevel verbosity) const
{
    using std::endl;

    const std::string tab = "   ";


    out << "CartesianGridElementAccessor<" << dim_ << "> info:" << endl;
    out.push(tab);

    if (contains(verbosity,VerbosityLevel::debug))
    	out << "Memory address = " << &(*this) << endl;

    CartesianGridElement<dim_>::print_info(out,verbosity);

    if (contains(verbosity,VerbosityLevel::debug))
    {
    	out << "Element cache memory address = " << &elem_values_ << endl;
    	elem_values_.print_info(out);

    	for (int i = 0 ; i < n_faces ; ++i)
    		out << "Face[" << i << "] cache memory address = " << &face_values_[i] << endl;
    }
    out.pop();

}


template <int dim_>
void
CartesianGridElementAccessor<dim_>::
ElementValuesCache::
print_info(LogStream &out) const
{
    using std::endl;

    const std::string tab = "   ";

    out << "ElementValuesCache info: (memory address = " << &(*this) << ")" << endl;
    out.push(tab);

    this->flags_handler_.print_info(out);

    out.pop();
}


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/cartesian_grid_element_accessor.inst>
